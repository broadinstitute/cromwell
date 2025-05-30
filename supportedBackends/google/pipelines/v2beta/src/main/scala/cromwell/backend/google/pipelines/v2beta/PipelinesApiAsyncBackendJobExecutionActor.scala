package cromwell.backend.google.pipelines.v2beta

import akka.actor.ActorRef
import cats.data.NonEmptyList
import cats.implicits._
import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import common.util.StringUtil._
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.io.PipelinesApiWorkingDisk
import cromwell.backend.google.pipelines.v2beta.PipelinesApiAsyncBackendJobExecutionActor._
import cromwell.backend.standard.StandardAsyncExecutionActorParams
import cromwell.backend.google.pipelines.common.PapiPollResultMonitorActor
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.filesystems.drs.DrsPath
import cromwell.filesystems.gcs.GcsPathBuilder.ValidFullGcsPath
import cromwell.filesystems.gcs.{GcsPath, GcsPathBuilder}
import org.apache.commons.codec.digest.DigestUtils
import org.apache.commons.csv.{CSVFormat, CSVPrinter}
import org.apache.commons.io.output.ByteArrayOutputStream
import wom.core.FullyQualifiedName
import wom.expression.FileEvaluation
import wom.values.{GlobFunctions, WomFile, WomGlobFile, WomSingleFile, WomUnlistedDirectory}

import java.nio.charset.Charset
import java.io.{FileNotFoundException, OutputStreamWriter}
import scala.concurrent.Future
import scala.io.Source
import scala.language.postfixOps
import scala.util.control.NoStackTrace

class PipelinesApiAsyncBackendJobExecutionActor(standardParams: StandardAsyncExecutionActorParams)
    extends cromwell.backend.google.pipelines.common.PipelinesApiAsyncBackendJobExecutionActor(standardParams)
    with PipelinesApiReferenceFilesMappingOperations {

  jobLogger.info("AN-436 Instantiated class")

  // The original implementation assumes the WomFiles are all WomMaybePopulatedFiles and wraps everything in a PipelinesApiFileInput
  // In v2 we can differentiate files from directories
  override protected def pipelinesApiInputsFromWomFiles(inputName: String,
                                                        remotePathArray: Seq[WomFile],
                                                        localPathArray: Seq[WomFile],
                                                        jobDescriptor: BackendJobDescriptor
  ): Iterable[PipelinesApiInput] =
    (remotePathArray zip localPathArray) flatMap {
      case (remotePath: WomUnlistedDirectory, localPath) =>
        Seq(
          PipelinesApiDirectoryInput(inputName,
                                     getPath(remotePath.valueString).get,
                                     DefaultPathBuilder.get(localPath.valueString),
                                     workingDisk
          )
        )
      case (remotePath, localPath) =>
        Seq(
          PipelinesApiFileInput(inputName,
                                getPath(remotePath.valueString).get,
                                DefaultPathBuilder.get(localPath.valueString),
                                workingDisk
          )
        )
    }

  // The original implementation recursively finds all non directory files, in V2 we can keep directory as is
  override protected lazy val callInputFiles: Map[FullyQualifiedName, Seq[WomFile]] = jobDescriptor.localInputs map {
    case (key, womFile) =>
      key -> womFile.collectAsSeq {
        case womFile: WomFile if !inputsToNotLocalize.contains(womFile) => womFile
      }
  }

  override val pollingResultMonitorActor: Option[ActorRef] = Option(
    context.actorOf(
      PapiPollResultMonitorActor.props(serviceRegistryActor,
                                       workflowDescriptor,
                                       jobDescriptor,
                                       validatedRuntimeAttributes,
                                       platform,
                                       jobLogger
      )
    )
  )

  private lazy val gcsTransferLibrary =
    Source.fromInputStream(Thread.currentThread.getContextClassLoader.getResourceAsStream("gcs_transfer.sh")).mkString

  private def gcsLocalizationTransferBundle[T <: PipelinesApiInput](
    gcsTransferConfiguration: GcsTransferConfiguration
  )(bucket: String, inputs: NonEmptyList[T]): String = {
    val project = inputs.head.cloudPath.asInstanceOf[GcsPath].projectId
    val maxAttempts = gcsTransferConfiguration.transferAttempts

    // Split files and directories out so files can possibly benefit from a `gsutil -m cp -I ...` optimization
    // on a per-container-parent-directory basis.
    val (files, directories) = inputs.toList partition { _.isInstanceOf[PipelinesApiFileInput] }

    // Files with different names between cloud and container are not eligible for bulk copying.
    val (filesWithSameNames, filesWithDifferentNames) = files partition { f =>
      f.cloudPath.asInstanceOf[GcsPath].nioPath.getFileName.toString == f.containerPath.getFileName.toString
    }

    val filesByContainerParentDirectory = filesWithSameNames.groupBy(_.containerPath.parent.toString)
    // Deduplicate any inputs since parallel localization can't deal with this.
    val uniqueFilesByContainerParentDirectory = filesByContainerParentDirectory map { case (p, fs) => p -> fs.toSet }

    val filesWithSameNamesTransferBundles: List[String] = uniqueFilesByContainerParentDirectory.toList map {
      case (containerParent, filesWithSameParent) =>
        val arrayIdentifier = s"files_to_localize_" + DigestUtils.md5Hex(bucket + containerParent)
        val entries = filesWithSameParent.map(_.cloudPath) mkString ("\"", "\"\n|  \"", "\"")

        s"""
           |# Localize files from source bucket '$bucket' to container parent directory '$containerParent'.
           |$arrayIdentifier=(
           |  "$project"   # project to use if requester pays
           |  "$maxAttempts" # max transfer attempts
           |  "${containerParent.ensureSlashed}" # container parent directory
           |  $entries
           |)
           |
           |localize_files "$${$arrayIdentifier[@]}"
       """.stripMargin
    }

    val filesWithDifferentNamesTransferBundles = filesWithDifferentNames map { f =>
      val arrayIdentifier =
        s"singleton_file_to_localize_" + DigestUtils.md5Hex(f.cloudPath.pathAsString + f.containerPath.pathAsString)
      s"""
         |# Localize singleton file '${f.cloudPath.pathAsString}' to '${f.containerPath.pathAsString}'.
         |$arrayIdentifier=(
         |  "$project"
         |  "$maxAttempts"
         |  "${f.cloudPath}"
         |  "${f.containerPath}"
         |)
         |
         |localize_singleton_file "$${$arrayIdentifier[@]}"
       """.stripMargin
    }

    // Only write a transfer bundle for directories if there are directories to be localized. Emptiness isn't a concern
    // for files since there is always at least the command script to be localized.
    val directoryTransferBundle =
      if (directories.isEmpty) ""
      else {
        val entries = directories flatMap { i => List(i.cloudPath, i.containerPath) } mkString ("\"", "\"\n|  \"", "\"")

        val arrayIdentifier = s"directories_to_localize_" + DigestUtils.md5Hex(bucket)

        s"""
           |# Directories from source bucket '$bucket'.
           |$arrayIdentifier=(
           |  "$project"    # project to use if requester pays
           |  "$maxAttempts" # max transfer attempts
           |  $entries
           |)
           |
           |localize_directories "$${$arrayIdentifier[@]}"
       """.stripMargin
      }

    (directoryTransferBundle :: (filesWithSameNamesTransferBundles ++ filesWithDifferentNamesTransferBundles)) mkString "\n\n"
  }

  private def gcsDelocalizationTransferBundle[T <: PipelinesApiOutput](
    transferConfiguration: GcsTransferConfiguration
  )(bucket: String, outputs: NonEmptyList[T]): String = {
    val project = outputs.head.cloudPath.asInstanceOf[GcsPath].projectId
    val maxAttempts = transferConfiguration.transferAttempts

    val transferItems = outputs.toList.flatMap { output =>
      val kind = output match {
        case o: PipelinesApiFileOutput if o.secondary => "file_or_directory" // if secondary the type is unknown
        case _: PipelinesApiFileOutput => "file" // a primary file
        case _: PipelinesApiDirectoryOutput => "directory" // a primary directory
      }

      val optional = Option(output) collectFirst {
        case o: PipelinesApiFileOutput if o.secondary || o.optional => "optional"
      } getOrElse "required"
      val contentType = output.contentType.map(_.toString).getOrElse("")

      List(kind, output.cloudPath.toString, output.containerPath.toString, optional, contentType)
    } mkString ("\"", "\"\n|  \"", "\"")

    val parallelCompositeUploadThreshold = jobDescriptor.workflowDescriptor.workflowOptions
      .getOrElse("parallel_composite_upload_threshold", transferConfiguration.parallelCompositeUploadThreshold)

    // Use a digest as bucket names can contain characters that are not legal in bash identifiers.
    val arrayIdentifier = s"delocalize_" + DigestUtils.md5Hex(bucket)
    s"""
       |# $bucket
       |$arrayIdentifier=(
       |  "$project"       # project
       |  "$maxAttempts"   # max attempts
       |  "$parallelCompositeUploadThreshold" # parallel composite upload threshold, will not be used for directory types
       |  $transferItems
       |)
       |
       |delocalize "$${$arrayIdentifier[@]}"
      """.stripMargin
  }

  private def bracketTransfersWithMessages(activity: String)(transferBody: String): String =
    List(
      s"timestamped_message '$activity script execution started...'",
      transferBody,
      s"timestamped_message '$activity script execution complete.'"
    ) mkString "\n"

  import mouse.all._

  override def uploadDrsLocalizationManifest(createPipelineParameters: CreatePipelineParameters,
                                             cloudPath: Path
  ): Future[Unit] = {
    val content = generateDrsLocalizerManifest(createPipelineParameters.inputOutputParameters.fileInputParameters)
    if (content.nonEmpty)
      asyncIo.writeAsync(cloudPath, content, Seq(CloudStorageOptions.withMimeType("text/plain")))
    else
      Future.unit
  }

  private def generateGcsLocalizationScript(inputs: List[PipelinesApiInput],
                                            referenceInputsToMountedPathsOpt: Option[Map[PipelinesApiInput, String]]
  )(implicit gcsTransferConfiguration: GcsTransferConfiguration): String = {
    // Generate a mapping of reference inputs to their mounted paths and a section of the localization script to
    // "faux localize" these reference inputs with symlinks to their locations on mounted reference disks.
    import cromwell.backend.google.pipelines.common.action.ActionUtils.shellEscaped
    val referenceFilesLocalizationScript = {
      val symlinkCreationCommandsOpt = referenceInputsToMountedPathsOpt map { referenceInputsToMountedPaths =>
        referenceInputsToMountedPaths map { case (input, absolutePathOnRefDisk) =>
          s"mkdir -p ${shellEscaped(input.containerPath.parent.pathAsString)} && ln -s ${shellEscaped(
              absolutePathOnRefDisk
            )} ${shellEscaped(input.containerPath.pathAsString)}"
        }
      }

      if (symlinkCreationCommandsOpt.exists(_.nonEmpty)) {
        s"""
           |# Faux-localizing reference files (if any) by creating symbolic links to the files located on the mounted reference disk
           |${symlinkCreationCommandsOpt.get.mkString("\n")}
           |""".stripMargin
      } else {
        "\n# No reference disks mounted / no symbolic links created since no matching reference files found in the inputs to this call.\n"
      }
    }

    val maybeReferenceFilesLocalizationScript =
      if (useReferenceDisks) {
        referenceFilesLocalizationScript
      } else {
        "\n# No reference disks mounted since not requested in workflow options.\n"
      }

    val regularFilesLocalizationScript = {
      val regularFiles = referenceInputsToMountedPathsOpt
        .map(maybeReferenceInputsToMountedPaths => inputs diff maybeReferenceInputsToMountedPaths.keySet.toList)
        .getOrElse(inputs)
      if (regularFiles.nonEmpty) {
        val bundleFunction = (gcsLocalizationTransferBundle(gcsTransferConfiguration) _).tupled
        generateGcsTransferScript(regularFiles, bundleFunction)
      } else {
        ""
      }
    }

    val combinedLocalizationScript =
      s"""
         |$maybeReferenceFilesLocalizationScript
         |
         |$regularFilesLocalizationScript
         |""".stripMargin

    combinedLocalizationScript |> bracketTransfersWithMessages("Localization")
  }

  private def generateGcsDelocalizationScript(
    outputs: List[PipelinesApiOutput]
  )(implicit gcsTransferConfiguration: GcsTransferConfiguration): String = {
    val bundleFunction = (gcsDelocalizationTransferBundle(gcsTransferConfiguration) _).tupled
    generateGcsTransferScript(outputs, bundleFunction) |> bracketTransfersWithMessages("Delocalization")
  }

  private def generateGcsTransferScript[T <: PipelinesParameter](items: List[T],
                                                                 bundleFunction: ((String, NonEmptyList[T])) => String
  ): String = {
    val gcsItems = items collect { case i if i.cloudPath.isInstanceOf[GcsPath] => i }
    groupParametersByGcsBucket(gcsItems) map bundleFunction mkString "\n"
  }

  override protected def uploadGcsTransferLibrary(createPipelineParameters: CreatePipelineParameters,
                                                  cloudPath: Path,
                                                  gcsTransferConfiguration: GcsTransferConfiguration
  ): Future[Unit] =
    asyncIo.writeAsync(cloudPath, gcsTransferLibrary, Seq(CloudStorageOptions.withMimeType("text/plain")))

  override def uploadGcsLocalizationScript(createPipelineParameters: CreatePipelineParameters,
                                           cloudPath: Path,
                                           transferLibraryContainerPath: Path,
                                           gcsTransferConfiguration: GcsTransferConfiguration,
                                           referenceInputsToMountedPathsOpt: Option[Map[PipelinesApiInput, String]]
  ): Future[Unit] = {
    val content = generateGcsLocalizationScript(createPipelineParameters.inputOutputParameters.fileInputParameters,
                                                referenceInputsToMountedPathsOpt
    )(gcsTransferConfiguration)
    asyncIo.writeAsync(cloudPath,
                       s"source '$transferLibraryContainerPath'\n\n" + content,
                       Seq(CloudStorageOptions.withMimeType("text/plain"))
    )
  }

  override def uploadGcsDelocalizationScript(createPipelineParameters: CreatePipelineParameters,
                                             cloudPath: Path,
                                             transferLibraryContainerPath: Path,
                                             gcsTransferConfiguration: GcsTransferConfiguration
  ): Future[Unit] = {
    val content = generateGcsDelocalizationScript(createPipelineParameters.inputOutputParameters.fileOutputParameters)(
      gcsTransferConfiguration
    )
    asyncIo.writeAsync(cloudPath,
                       s"source '$transferLibraryContainerPath'\n\n" + content,
                       Seq(CloudStorageOptions.withMimeType("text/plain"))
    )
  }

  // Simply create a PipelinesApiDirectoryOutput in v2 instead of globbing
  override protected def generateUnlistedDirectoryOutputs(unlistedDirectory: WomUnlistedDirectory,
                                                          fileEvaluation: FileEvaluation
  ): List[PipelinesApiOutput] = {
    val destination = callRootPath.resolve(unlistedDirectory.value.stripPrefix("/"))
    val (relpath, disk) = relativePathAndAttachedDisk(unlistedDirectory.value, runtimeAttributes.disks)
    val directoryOutput = PipelinesApiDirectoryOutput(makeSafeReferenceName(unlistedDirectory.value),
                                                      destination,
                                                      relpath,
                                                      disk,
                                                      fileEvaluation.optional,
                                                      fileEvaluation.secondary
    )
    List(directoryOutput)
  }

  // De-localize the glob directory as a PipelinesApiDirectoryOutput instead of using * pattern match
  override def generateGlobFileOutputs(womFile: WomGlobFile): List[PipelinesApiOutput] = {
    val globName = GlobFunctions.globName(womFile.value)
    val globDirectory = globName + "/"
    val globListFile = globName + ".list"
    val gcsGlobDirectoryDestinationPath = callRootPath.resolve(globDirectory)
    val gcsGlobListFileDestinationPath = callRootPath.resolve(globListFile)

    val (_, globDirectoryDisk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)

    // We need both the glob directory and the glob list:
    List(
      // The glob directory:
      PipelinesApiDirectoryOutput(
        makeSafeReferenceName(globDirectory),
        gcsGlobDirectoryDestinationPath,
        DefaultPathBuilder.get(globDirectory),
        globDirectoryDisk,
        optional = false,
        secondary = false
      ),
      // The glob list file:
      PipelinesApiFileOutput(
        makeSafeReferenceName(globListFile),
        gcsGlobListFileDestinationPath,
        DefaultPathBuilder.get(globListFile),
        globDirectoryDisk,
        optional = false,
        secondary = false
      )
    )
  }

  override def womFileToGcsPath(jesOutputs: Set[PipelinesApiOutput])(womFile: WomFile): WomFile =
    womFile mapFile { path =>
      jesOutputs collectFirst {
        case jesOutput if jesOutput.name == makeSafeReferenceName(path) =>
          val pathAsString = jesOutput.cloudPath.pathAsString
          if (jesOutput.isFileParameter && !jesOutput.cloudPath.exists) {
            // This is not an error if the path represents a `File?` optional output (the PAPI delocalization script
            // should have failed if this file output was not optional but missing). Throw to produce the correct "empty
            // optional" value for a missing optional file output.
            throw new FileNotFoundException(s"GCS output file not found: $pathAsString")
          }
          pathAsString
      } getOrElse {
        GcsPathBuilder.validateGcsPath(path) match {
          case _: ValidFullGcsPath => path

          /*
           * Strip the prefixes in RuntimeOutputMapping.prefixFilters from the path, one at a time.
           * For instance
           * file:///cromwell_root/bucket/workflow_name/6d777414-5ee7-4c60-8b9e-a02ec44c398e/call-A/file.txt will progressively become
           *
           * /cromwell_root/bucket/workflow_name/6d777414-5ee7-4c60-8b9e-a02ec44c398e/call-A/file.txt
           * bucket/workflow_name/6d777414-5ee7-4c60-8b9e-a02ec44c398e/call-A/file.txt
           * call-A/file.txt
           *
           * This code is called as part of a path mapper that will be applied to the WOMified cwl.output.json.
           * The cwl.output.json when it's being read by Cromwell from the bucket still contains local paths
           * (as they were created by the cwl tool).
           * In order to keep things working we need to map those local paths to where they were actually delocalized,
           * which is determined in cromwell.backend.google.pipelines.v2beta.api.Delocalization.
           */
          case _ =>
            (callRootPath /
              RuntimeOutputMapping
                .prefixFilters(workflowPaths.workflowRoot)
                .foldLeft(path) { case (newPath, prefix) =>
                  newPath.stripPrefix(prefix)
                }).pathAsString
        }
      }
    }

  override def generateSingleFileOutputs(womFile: WomSingleFile,
                                         fileEvaluation: FileEvaluation
  ): List[PipelinesApiFileOutput] = {
    val (relpath, disk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)
    // If the file is on a custom mount point, resolve it so that the full mount path will show up in the cloud path
    // For the default one (cromwell_root), the expectation is that it does not appear
    val mountedPath =
      if (!disk.mountPoint.isSamePathAs(PipelinesApiWorkingDisk.Default.mountPoint)) disk.mountPoint.resolve(relpath)
      else relpath
    // Normalize the local path (to get rid of ".." and "."). Also strip any potential leading / so that it gets appended to the call root
    val normalizedPath = mountedPath.normalize().pathAsString.stripPrefix("/")
    val destination = callRootPath.resolve(normalizedPath)
    val jesFileOutput = PipelinesApiFileOutput(makeSafeReferenceName(womFile.value),
                                               destination,
                                               relpath,
                                               disk,
                                               fileEvaluation.optional,
                                               fileEvaluation.secondary
    )
    List(jesFileOutput)
  }
}

object PipelinesApiAsyncBackendJobExecutionActor {
  // GCS path regexes comments:
  // - The (?s) option at the start makes '.' expression to match any symbol, including '\n'
  // - All valid GCS paths start with gs://
  // - Bucket names:
  //  - The bucket name is matched inside a set of '()'s so it can be used later.
  //  - The bucket name must start with a letter or number (https://cloud.google.com/storage/docs/naming)
  //  - Then, anything that is not a '/' is part of the bucket name
  // - Allow zero or more subdirectories, with (/[^/]+)*
  // - Then, for files:
  //  - There must be at least one '/', followed by some content in the file name.
  // - Or, then, for directories:
  //  - If we got this far, we already have a valid directory path. Allow it to optionally end with a `/` character.
  private val gcsFilePathMatcher = "(?s)^gs://([a-zA-Z0-9][^/]+)(/[^/]+)*/[^/]+$".r
  private val gcsDirectoryPathMatcher = "(?s)^gs://([a-zA-Z0-9][^/]+)(/[^/]+)*/?$".r

  private[v2beta] def groupParametersByGcsBucket[T <: PipelinesParameter](
    parameters: List[T]
  ): Map[String, NonEmptyList[T]] =
    parameters.map { param =>
      def pathTypeString = if (param.isFileParameter) "File" else "Directory"
      val regexToUse = if (param.isFileParameter) gcsFilePathMatcher else gcsDirectoryPathMatcher

      param.cloudPath.pathAsString match {
        case regexToUse(bucket) => Map(bucket -> NonEmptyList.of(param))
        case regexToUse(bucket, _) => Map(bucket -> NonEmptyList.of(param))
        case other =>
          throw new Exception(
            s"$pathTypeString path '$other' did not match the expected regex: ${regexToUse.pattern.toString}"
          ) with NoStackTrace
      }
    } combineAll

  private[v2beta] def generateDrsLocalizerManifest(inputs: List[PipelinesApiInput]): String = {
    val outputStream = new ByteArrayOutputStream()
    val csvPrinter = new CSVPrinter(new OutputStreamWriter(outputStream), CSVFormat.DEFAULT)
    val drsFileInputs = inputs collect { case drsInput @ PipelinesApiFileInput(_, drsPath: DrsPath, _, _) =>
      (drsInput, drsPath)
    }
    drsFileInputs foreach { case (drsInput, drsPath) =>
      csvPrinter.printRecord(drsPath.pathAsString, drsInput.containerPath.pathAsString)
    }
    csvPrinter.close(true)
    outputStream.toString(Charset.defaultCharset())
  }
}
