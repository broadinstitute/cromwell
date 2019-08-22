package cromwell.backend.google.pipelines.v2alpha1

import cats.data.NonEmptyList
import cats.instances.list._
import cats.instances.map._
import cats.syntax.foldable._
import com.google.api.services.genomics.v2alpha1.model.Mount
import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import cromwell.backend.BackendJobDescriptor
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.common.io.PipelinesApiWorkingDisk
import cromwell.backend.google.pipelines.v2alpha1.PipelinesApiAsyncBackendJobExecutionActor._
import cromwell.backend.standard.StandardAsyncExecutionActorParams
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.filesystems.gcs.GcsPathBuilder.ValidFullGcsPath
import cromwell.filesystems.gcs.{GcsPath, GcsPathBuilder}
import org.apache.commons.codec.digest.DigestUtils
import wom.core.FullyQualifiedName
import wom.expression.FileEvaluation
import wom.values.{GlobFunctions, WomFile, WomGlobFile, WomMaybeListedDirectory, WomMaybePopulatedFile, WomSingleFile, WomUnlistedDirectory}

import scala.concurrent.Future
import scala.io.Source
import scala.language.postfixOps

class PipelinesApiAsyncBackendJobExecutionActor(standardParams: StandardAsyncExecutionActorParams) extends cromwell.backend.google.pipelines.common.PipelinesApiAsyncBackendJobExecutionActor(standardParams) {

  // The original implementation assumes the WomFiles are all WomMaybePopulatedFiles and wraps everything in a PipelinesApiFileInput
  // In v2 we can differentiate files from directories 
  override protected def pipelinesApiInputsFromWomFiles(inputName: String,
                                                        remotePathArray: Seq[WomFile],
                                                        localPathArray: Seq[WomFile],
                                                        jobDescriptor: BackendJobDescriptor): Iterable[PipelinesApiInput] = {
    (remotePathArray zip localPathArray) flatMap {
      case (remotePath: WomMaybeListedDirectory, localPath) =>
        maybeListedDirectoryToPipelinesParameters(inputName, remotePath, localPath.valueString)
      case (remotePath: WomUnlistedDirectory, localPath) =>
        Seq(PipelinesApiDirectoryInput(inputName, getPath(remotePath.valueString).get, DefaultPathBuilder.get(localPath.valueString), workingDisk))
      case (remotePath: WomMaybePopulatedFile, localPath) =>
        maybePopulatedFileToPipelinesParameters(inputName, remotePath, localPath.valueString)
      case (remotePath, localPath) =>
        Seq(PipelinesApiFileInput(inputName, getPath(remotePath.valueString).get, DefaultPathBuilder.get(localPath.valueString), workingDisk))
    }
  }

  // The original implementation recursively finds all non directory files, in V2 we can keep directory as is
  override protected lazy val callInputFiles: Map[FullyQualifiedName, Seq[WomFile]] = jobDescriptor.localInputs map {
    case (key, womFile) =>
      key -> womFile.collectAsSeq({
        case womFile: WomFile if !inputsToNotLocalize.contains(womFile) => womFile
      })
  }

  private lazy val transferScriptTemplate =
    Source.fromInputStream(Thread.currentThread.getContextClassLoader.getResourceAsStream("gcs_transfer.sh")).mkString

  private def gcsLocalizationTransferBundle[T <: PipelinesApiInput](localizationConfiguration: LocalizationConfiguration)(bucket: String, inputs: NonEmptyList[T]): String = {
    val project = inputs.head.cloudPath.asInstanceOf[GcsPath].projectId
    val maxAttempts = localizationConfiguration.localizationAttempts
    val transferItems = inputs.toList.flatMap { i =>
      val kind = i match {
        case _: PipelinesApiFileInput => "file"
        case _: PipelinesApiDirectoryInput => "directory"
      }
      List(kind, i.cloudPath, i.containerPath)
    } mkString("\"", "\"\n|  \"", "\"")

    // Use a digest as bucket names can contain characters that are not legal in bash identifiers.
    val arrayIdentifier = s"localize_" + DigestUtils.md5Hex(bucket)
    s"""
       |# $bucket
       |$arrayIdentifier=(
       |  "localize" # direction
       |  "$project"       # project
       |  "$maxAttempts"   # max attempts
       |  $transferItems
       |)
       |
       |transfer "$${$arrayIdentifier[@]}"
      """.stripMargin
  }

  private def gcsDelocalizationTransferBundle[T <: PipelinesApiOutput](localizationConfiguration: LocalizationConfiguration)(bucket: String, outputs: NonEmptyList[T]): String = {
    val project = outputs.head.cloudPath.asInstanceOf[GcsPath].projectId
    val maxAttempts = localizationConfiguration.localizationAttempts

    val transferItems = outputs.toList.flatMap { output =>
      val kind = output match {
        case o: PipelinesApiFileOutput if o.secondary => "file_or_directory" // if secondary the type is unknown
        case _: PipelinesApiFileOutput => "file" // a primary file
        case _: PipelinesApiDirectoryOutput => "directory" // a primary directory
      }

      val optional = Option(output) collectFirst { case o: PipelinesApiFileOutput if o.secondary || o.optional => "optional" } getOrElse "required"
      val contentType = output.contentType.getOrElse("")

      List(kind, output.cloudPath, output.containerPath, optional, contentType)
    } mkString("\"", "\"\n|  \"", "\"")


    // Use a digest as bucket names can contain characters that are not legal in bash identifiers.
    val arrayIdentifier = s"delocalize_" + DigestUtils.md5Hex(bucket)
    s"""
       |# $bucket
       |$arrayIdentifier=(
       |  "delocalize" # direction
       |  "$project"       # project
       |  "$maxAttempts"   # max attempts
       |  $transferItems
       |)
       |
       |transfer "$${$arrayIdentifier[@]}"
      """.stripMargin
  }

  private def generateGcsLocalizationScript(inputs: List[PipelinesApiInput], mounts: List[Mount])(implicit localizationConfiguration: LocalizationConfiguration): String = {
    val bundleFunction = (gcsLocalizationTransferBundle(localizationConfiguration) _).tupled
    generateGcsTransferScript(inputs, mounts, bundleFunction)
  }

  private def generateGcsDelocalizationScript(outputs: List[PipelinesApiOutput], mounts: List[Mount])(implicit localizationConfiguration: LocalizationConfiguration): String = {
    val bundleFunction = (gcsDelocalizationTransferBundle(localizationConfiguration) _).tupled
    generateGcsTransferScript(outputs, mounts, bundleFunction)
  }

  private def generateGcsTransferScript[T <: PipelinesParameter](items: List[T], mounts: List[Mount], bundleFunction: ((String, NonEmptyList[T])) => String): String = {
    val gcsItems = items collect { case i if i.cloudPath.isInstanceOf[GcsPath] => i }
    val localizationBundles = groupParametersByGcsBucket(gcsItems) map bundleFunction mkString "\n"
    transferScriptTemplate + localizationBundles
  }

  override def uploadGcsLocalizationScript(createPipelineParameters: CreatePipelineParameters, cloudPath: Path, localizationConfiguration: LocalizationConfiguration): Future[Unit] = {
    val mounts = PipelinesConversions.toMounts(createPipelineParameters)
    val content = generateGcsLocalizationScript(createPipelineParameters.inputOutputParameters.fileInputParameters, mounts)(localizationConfiguration)
    asyncIo.writeAsync(cloudPath, content, Seq(CloudStorageOptions.withMimeType("text/plain")))
  }

  override def uploadGcsDelocalizationScript(createPipelineParameters: CreatePipelineParameters, cloudPath: Path, localizationConfiguration: LocalizationConfiguration): Future[Unit] = {
    val mounts = PipelinesConversions.toMounts(createPipelineParameters)
    val content = generateGcsDelocalizationScript(createPipelineParameters.inputOutputParameters.fileOutputParameters, mounts)(localizationConfiguration)
    asyncIo.writeAsync(cloudPath, content, Seq(CloudStorageOptions.withMimeType("text/plain")))
  }

  // Simply create a PipelinesApiDirectoryOutput in v2 instead of globbing
  override protected def generateUnlistedDirectoryOutputs(unlistedDirectory: WomUnlistedDirectory, fileEvaluation: FileEvaluation): List[PipelinesApiOutput] = {
    val destination = callRootPath.resolve(unlistedDirectory.value.stripPrefix("/"))
    val (relpath, disk) = relativePathAndAttachedDisk(unlistedDirectory.value, runtimeAttributes.disks)
    val directoryOutput = PipelinesApiDirectoryOutput(makeSafeReferenceName(unlistedDirectory.value), destination, relpath, disk, fileEvaluation.optional, fileEvaluation.secondary)
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
      PipelinesApiDirectoryOutput(makeSafeReferenceName(globDirectory), gcsGlobDirectoryDestinationPath, DefaultPathBuilder.get(globDirectory), globDirectoryDisk, optional = false, secondary = false),
      // The glob list file:
      PipelinesApiFileOutput(makeSafeReferenceName(globListFile), gcsGlobListFileDestinationPath, DefaultPathBuilder.get(globListFile), globDirectoryDisk, optional = false, secondary = false)
    )
  }

  override def womFileToGcsPath(jesOutputs: Set[PipelinesApiOutput])(womFile: WomFile): WomFile = {
    womFile mapFile { path =>
      jesOutputs collectFirst {
        case jesOutput if jesOutput.name == makeSafeReferenceName(path) => jesOutput.cloudPath.pathAsString
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
            * which is determined in cromwell.backend.google.pipelines.v2alpha1.api.Delocalization.
            */
          case _ => (callRootPath /
            RuntimeOutputMapping
                .prefixFilters(workflowPaths.workflowRoot)
                .foldLeft(path)({
                  case (newPath, prefix) => newPath.stripPrefix(prefix)
                })
            ).pathAsString
        }
      }
    }
  }

  private def maybePopulatedFileToPipelinesParameters(inputName: String, maybePopulatedFile: WomMaybePopulatedFile, localPath: String) = {
    val secondaryFiles = maybePopulatedFile.secondaryFiles.flatMap({ secondaryFile =>
      pipelinesApiInputsFromWomFiles(secondaryFile.valueString, List(secondaryFile), List(relativeLocalizationPath(secondaryFile)), jobDescriptor)
    })

    Seq(PipelinesApiFileInput(inputName, getPath(maybePopulatedFile.valueString).get, DefaultPathBuilder.get(localPath), workingDisk)) ++ secondaryFiles
  }

  private def maybeListedDirectoryToPipelinesParameters(inputName: String, womMaybeListedDirectory: WomMaybeListedDirectory, localPath: String) = womMaybeListedDirectory match {
    // If there is a path, simply localize as a directory
    case WomMaybeListedDirectory(Some(path), _, _, _) =>
      List(PipelinesApiDirectoryInput(inputName, getPath(path).get, DefaultPathBuilder.get(localPath), workingDisk))

    // If there is a listing, recurse and call pipelinesApiInputsFromWomFiles on all the listed files
    case WomMaybeListedDirectory(_, Some(listing), _, _) if listing.nonEmpty =>
      listing.flatMap({
        case womFile: WomFile if isAdHocFile(womFile) =>
          pipelinesApiInputsFromWomFiles(makeSafeReferenceName(womFile.valueString), List(womFile), List(fileName(womFile)), jobDescriptor)
        case womFile: WomFile =>
          pipelinesApiInputsFromWomFiles(makeSafeReferenceName(womFile.valueString), List(womFile), List(relativeLocalizationPath(womFile)), jobDescriptor)
      })
    case _ => List.empty
  }

  override def generateSingleFileOutputs(womFile: WomSingleFile, fileEvaluation: FileEvaluation) = {
    val (relpath, disk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)
    // If the file is on a custom mount point, resolve it so that the full mount path will show up in the cloud path
    // For the default one (cromwell_root), the expectation is that it does not appear
    val mountedPath = if (!disk.mountPoint.isSamePathAs(PipelinesApiWorkingDisk.Default.mountPoint)) disk.mountPoint.resolve(relpath) else relpath
    // Normalize the local path (to get rid of ".." and "."). Also strip any potential leading / so that it gets appended to the call root
    val normalizedPath = mountedPath.normalize().pathAsString.stripPrefix("/")
    val destination = callRootPath.resolve(normalizedPath)
    val jesFileOutput = PipelinesApiFileOutput(makeSafeReferenceName(womFile.value), destination, relpath, disk, fileEvaluation.optional, fileEvaluation.secondary)
    List(jesFileOutput)
  }
}

object PipelinesApiAsyncBackendJobExecutionActor {
  private val gcsPathMatcher = "^gs://?([^/]+)/.*".r

  private [v2alpha1] def groupParametersByGcsBucket[T <: PipelinesParameter](parameters: List[T]): Map[String, NonEmptyList[T]] = {
    parameters.map { param =>
      param.cloudPath.toString match {
        case gcsPathMatcher(bucket) =>
          Map(bucket -> NonEmptyList.of(param))
        case other =>
          throw new RuntimeException(s"Programmer error: Path '$other' did not match the expected regex. This should have been impossible")
      }
    } combineAll
  }
}
