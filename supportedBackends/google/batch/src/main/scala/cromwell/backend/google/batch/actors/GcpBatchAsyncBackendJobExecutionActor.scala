package cromwell.backend.google.batch.actors

import _root_.io.grpc.{Status => GrpcStatus}
import akka.actor.{ActorLogging, ActorRef}
import akka.http.scaladsl.model.{ContentType, ContentTypes}
import akka.pattern.AskSupport
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.implicits._
import com.google.cloud.batch.v1.JobName
import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import common.util.StringUtil._
import common.validation.ErrorOr.ErrorOr
import cromwell.backend._
import cromwell.backend.async.{
  AbortedExecutionHandle,
  ExecutionHandle,
  FailedNonRetryableExecutionHandle,
  FailedRetryableExecutionHandle,
  PendingExecutionHandle
}
import cromwell.backend.google.batch.GcpBatchBackendLifecycleActorFactory
import cromwell.backend.google.batch.actors.BatchApiRunCreationClient.JobAbortedException
import cromwell.backend.google.batch.api.GcpBatchRequestFactory._
import cromwell.backend.google.batch.errors.FailedToDelocalizeFailure
import cromwell.backend.google.batch.io._
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.batch.models.GcpBatchJobPaths.GcsTransferLibraryName
import cromwell.backend.google.batch.models.RunStatus.TerminalRunStatus
import cromwell.backend.google.batch.models._
import cromwell.backend.google.batch.monitoring.{BatchInstrumentation, CheckpointingConfiguration, MonitoringImage}
import cromwell.backend.google.batch.runnable.WorkflowOptionKeys
import cromwell.backend.google.batch.util.{GcpBatchReferenceFilesMappingOperations, RuntimeOutputMapping}
import cromwell.filesystems.gcs.GcsPathBuilder
import cromwell.filesystems.gcs.GcsPathBuilder.ValidFullGcsPath

import java.io.FileNotFoundException
import cromwell.backend.standard.{
  ScriptPreambleData,
  StandardAdHocValue,
  StandardAsyncExecutionActor,
  StandardAsyncExecutionActorParams,
  StandardAsyncJob
}
import cromwell.core._
import cromwell.core.io.IoCommandBuilder
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.filesystems.drs.{DrsPath, DrsResolver}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.http.HttpPath
import cromwell.filesystems.sra.SraPath
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.keyvalue.KeyValueServiceActor.{KvJobKey, KvPair, ScopedKey}
import cromwell.services.metadata.CallMetadataKeys
import mouse.all._
import shapeless.Coproduct
import org.apache.commons.codec.digest.DigestUtils
import org.apache.commons.csv.{CSVFormat, CSVPrinter}
import org.apache.commons.io.output.ByteArrayOutputStream
import wom.callable.Callable.OutputDefinition
import wom.callable.MetaValueElement.{MetaValueElementBoolean, MetaValueElementObject}
import wom.callable.AdHocValue
import wom.core.FullyQualifiedName
import wom.expression.{FileEvaluation, NoIoFunctionSet}
import wom.values._

import java.io.OutputStreamWriter
import java.net.SocketTimeoutException
import java.nio.charset.Charset
import java.util.Base64
import scala.concurrent.Future
import scala.concurrent.duration._
import scala.io.Source
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scala.util.control.NoStackTrace

object GcpBatchAsyncBackendJobExecutionActor {
  // TODO: Clean these up
  val maxUnexpectedRetries = 2

  val JesFailedToDelocalize = 5
  val JesUnexpectedTermination = 13
  val JesPreemption = 14

  val BatchFailedPreConditionErrorCode = 9
  val BatchMysteriouslyCrashedErrorCode = 10

  // If the JES code is 2 (UNKNOWN), this sub-string indicates preemption:
  val FailedToStartDueToPreemptionSubstring = "failed to start due to preemption"
  val FailedV2Style = "The assigned worker has failed to complete the operation"

  def StandardException(errorCode: GrpcStatus,
                        message: String,
                        jobTag: String,
                        returnCodeOption: Option[Int],
                        stderrPath: Path
  ): Exception = {
    val returnCodeMessage = returnCodeOption match {
      case Some(returnCode) if returnCode == 0 => "Job exited without an error, exit code 0."
      case Some(returnCode) => s"Job exit code $returnCode. Check $stderrPath for more information."
      case None => "The job was stopped before the command finished."
    }

    new Exception(s"Task $jobTag failed. $returnCodeMessage Batch error code ${errorCode.getCode.value}. $message")
  }

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

  val GcpBatchOperationIdKey = "__gcp_batch_operation_id"

  type GcpBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, Run, RunStatus]

  val plainTextContentType: Option[ContentType.WithCharset] = Option(ContentTypes.`text/plain(UTF-8)`)

  private[batch] def groupParametersByGcsBucket[T <: BatchParameter](
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

  private[batch] def generateDrsLocalizerManifest(inputs: List[GcpBatchInput]): String = {
    val outputStream = new ByteArrayOutputStream()
    val csvPrinter = new CSVPrinter(new OutputStreamWriter(outputStream), CSVFormat.DEFAULT)
    val drsFileInputs = inputs collect { case drsInput @ GcpBatchFileInput(_, drsPath: DrsPath, _, _) =>
      (drsInput, drsPath)
    }
    drsFileInputs foreach { case (drsInput, drsPath) =>
      csvPrinter.printRecord(drsPath.pathAsString, drsInput.containerPath.pathAsString)
    }
    csvPrinter.close(true)
    outputStream.toString(Charset.defaultCharset())
  }

}

class GcpBatchAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
    extends BackendJobLifecycleActor
    with StandardAsyncExecutionActor
    with BatchApiRunCreationClient
    with BatchApiFetchJobClient
    with BatchApiAbortClient
    with AskSupport
    with GcpBatchJobCachingActorHelper
    with GcpBatchReferenceFilesMappingOperations
    with BatchInstrumentation
    with ActorLogging
    with CromwellInstrumentation {

  import GcpBatchAsyncBackendJobExecutionActor._
  override lazy val ioCommandBuilder: IoCommandBuilder = GcsBatchCommandBuilder

  lazy val workflowId: WorkflowId = jobDescriptor.workflowDescriptor.id

  /** The type of the run info when a job is started. */
  override type StandardAsyncRunInfo = Run

  /** The type of the run status returned during each poll. */
  override type StandardAsyncRunState = RunStatus

  override def receive: Receive =
    runCreationClientReceive orElse pollingActorClientReceive orElse abortActorClientReceive orElse kvClientReceive orElse super.receive

  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
   * in the state type.
   */
  def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean = thiz == that

  protected lazy val cmdInput: GcpBatchFileInput =
    GcpBatchFileInput(GcpBatchJobPaths.BatchExecParamName,
                      gcpBatchCallPaths.script,
                      DefaultPathBuilder.get(gcpBatchCallPaths.scriptFilename),
                      workingDisk
    )

  private lazy val jobDockerImage = jobDescriptor.maybeCallCachingEligible.dockerHash
    .getOrElse(runtimeAttributes.dockerImage)

  override def dockerImageUsed: Option[String] = Option(jobDockerImage)

  // Need to add previousRetryReasons and preemptible in order to get preemptible to work in the tests
  protected val previousRetryReasons: ErrorOr[PreviousRetryReasons] =
    PreviousRetryReasons.tryApply(jobDescriptor.prefetchedKvStoreEntries, jobDescriptor.key.attempt)

  override lazy val preemptible: Boolean = previousRetryReasons match {
    case Valid(PreviousRetryReasons(p, _)) => p < maxPreemption
    case _ => false
  }

  override def tryAbort(job: StandardAsyncJob): Unit =
    abortJob(workflowId = workflowId,
             jobName = JobName.parse(job.jobId),
             backendSingletonActor = backendSingletonActor,
             requestFactory = initializationData.requestFactory
    )

  val backendSingletonActor: ActorRef = standardParams.backendSingletonActorOption
    .getOrElse(throw new RuntimeException("GCP Batch actor cannot exist without its backend singleton 2"))

  /**
   * Takes two arrays of remote and local WOM File paths and generates the necessary `GcpBatchInput`s.
   */
  protected def gcpBatchInputsFromWomFiles(inputName: String,
                                           remotePathArray: Seq[WomFile],
                                           localPathArray: Seq[WomFile],
                                           jobDescriptor: BackendJobDescriptor
  ): Iterable[GcpBatchInput] =
    (remotePathArray zip localPathArray) flatMap {
      case (remotePath: WomMaybeListedDirectory, localPath) =>
        maybeListedDirectoryToBatchParameters(inputName, remotePath, localPath.valueString)
      case (remotePath: WomUnlistedDirectory, localPath) =>
        Seq(
          GcpBatchDirectoryInput(inputName,
                                 getPath(remotePath.valueString).get,
                                 DefaultPathBuilder.get(localPath.valueString),
                                 workingDisk
          )
        )
      case (remotePath: WomMaybePopulatedFile, localPath) =>
        maybePopulatedFileToBatchParameters(inputName, remotePath, localPath.valueString)
      case (remotePath, localPath) =>
        Seq(
          GcpBatchFileInput(inputName,
                            getPath(remotePath.valueString).get,
                            DefaultPathBuilder.get(localPath.valueString),
                            workingDisk
          )
        )
    }

  private def maybePopulatedFileToBatchParameters(inputName: String,
                                                  maybePopulatedFile: WomMaybePopulatedFile,
                                                  localPath: String
  ) = {
    val secondaryFiles = maybePopulatedFile.secondaryFiles.flatMap { secondaryFile =>
      gcpBatchInputsFromWomFiles(secondaryFile.valueString,
                                 List(secondaryFile),
                                 List(relativeLocalizationPath(secondaryFile)),
                                 jobDescriptor
      )
    }

    Seq(
      GcpBatchFileInput(inputName,
                        getPath(maybePopulatedFile.valueString).get,
                        DefaultPathBuilder.get(localPath),
                        workingDisk
      )
    ) ++ secondaryFiles
  }

  /**
   * Turns WomFiles into relative paths.  These paths are relative to the working disk.
   *
   * relativeLocalizationPath("foo/bar.txt") -> "foo/bar.txt"
   * relativeLocalizationPath("gs://some/bucket/foo.txt") -> "some/bucket/foo.txt"
   */
  override protected def relativeLocalizationPath(file: WomFile): WomFile =
    file.mapFile(value =>
      getPath(value) match {
        case Success(drsPath: DrsPath) => DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
        case Success(path) => path.pathWithoutScheme
        case _ => value
      }
    )

  lazy val localMonitoringImageScriptPath: Path =
    DefaultPathBuilder.get(gcpBatchCallPaths.batchMonitoringImageScriptFilename)

  override protected def fileName(file: WomFile): WomFile =
    file.mapFile(value =>
      getPath(value) match {
        case Success(drsPath: DrsPath) =>
          DefaultPathBuilder
            .get(DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync())
            .name
        case Success(path) => path.name
        case _ => value
      }
    )

  override lazy val inputsToNotLocalize: Set[WomFile] = {
    val localizeOptional = jobDescriptor.findInputFilesByParameterMeta {
      case MetaValueElementObject(values) => values.get("localization_optional").contains(MetaValueElementBoolean(true))
      case _ => false
    }
    val localizeSkipped = localizeOptional.filter(canSkipLocalize)
    val localizeMapped = localizeSkipped.map(cloudResolveWomFile)
    localizeSkipped ++ localizeMapped
  }

  private def canSkipLocalize(womFile: WomFile): Boolean = {
    var canSkipLocalize = true
    womFile.mapFile { value =>
      getPath(value) match {
        case Success(drsPath: DrsPath) =>
          val gsUriOption = DrsResolver.getSimpleGsUri(drsPath).unsafeRunSync()
          if (gsUriOption.isEmpty) {
            canSkipLocalize = false
          }
        case _ => /* ignore */
      }
      value
    }
    canSkipLocalize
  }

  // The original implementation recursively finds all non directory files, in V2 we can keep directory as is
  protected lazy val callInputFiles: Map[FullyQualifiedName, Seq[WomFile]] =
    // NOTE: This causes the tests to fail
    jobDescriptor.localInputs map { case (key, womFile) =>
      key -> womFile.collectAsSeq {
        case womFile: WomFile if !inputsToNotLocalize.contains(womFile) => womFile
      }
    }

  private lazy val gcsTransferLibrary =
    Source.fromInputStream(Thread.currentThread.getContextClassLoader.getResourceAsStream("gcs_transfer.sh")).mkString

  private def gcsLocalizationTransferBundle[T <: GcpBatchInput](
    gcsTransferConfiguration: GcsTransferConfiguration
  )(bucket: String, inputs: NonEmptyList[T]): String = {
    val project = inputs.head.cloudPath.asInstanceOf[GcsPath].projectId
    val maxAttempts = gcsTransferConfiguration.transferAttempts

    // Split files and directories out so files can possibly benefit from a `gsutil -m cp -I ...` optimization
    // on a per-container-parent-directory basis.
    val (files, directories) = inputs.toList partition {
      _.isInstanceOf[GcpBatchFileInput]
    }

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

  private def gcsDelocalizationTransferBundle[T <: GcpBatchOutput](
    transferConfiguration: GcsTransferConfiguration
  )(bucket: String, outputs: NonEmptyList[T]): String = {
    val project = outputs.head.cloudPath.asInstanceOf[GcsPath].projectId
    val maxAttempts = transferConfiguration.transferAttempts

    val transferItems = outputs.toList.flatMap { output =>
      val kind = output match {
        case o: GcpBatchFileOutput if o.secondary => "file_or_directory" // if secondary the type is unknown
        case _: GcpBatchFileOutput => "file" // a primary file
        case _: GcpBatchDirectoryOutput => "directory" // a primary directory
      }

      val optional = Option(output) collectFirst {
        case o: GcpBatchFileOutput if o.secondary || o.optional => "optional"
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

  def uploadDrsLocalizationManifest(createParameters: CreateBatchJobParameters, cloudPath: Path): Future[Unit] = {
    val content = generateDrsLocalizerManifest(createParameters.inputOutputParameters.fileInputParameters)
    if (content.nonEmpty)
      asyncIo.writeAsync(cloudPath, content, Seq(CloudStorageOptions.withMimeType("text/plain")))
    else
      Future.unit
  }

  def uploadScriptFile(): Future[Unit] =
    commandScriptContents
      .fold(
        errors =>
          Future
            .failed(
              new RuntimeException(
                errors.toList
                  .mkString(", ")
              )
            ),
        asyncIo
          .writeAsync(jobPaths.script, _, Seq(CloudStorageOptions.withMimeType("text/plain")))
      )

  def sendGoogleLabelsToMetadata(customLabels: Seq[GcpLabel]): Unit = {
    lazy val backendLabelEvents: Map[String, String] =
      ((backendLabels ++ customLabels) map { l => s"${CallMetadataKeys.BackendLabels}:${l.key}" -> l.value }).toMap
    tellMetadata(backendLabelEvents)
  }

  protected val useReferenceDisks: Boolean = {
    val optionName = WorkflowOptions.UseReferenceDisks.name
    workflowDescriptor.workflowOptions.getBoolean(optionName) match {
      case Success(value) => value
      case Failure(OptionNotFoundException(_)) => false
      case Failure(f) =>
        // Should not happen, this case should have been screened for and fast-failed during workflow materialization.
        log.error(
          f,
          s"Programmer error: unexpected failure attempting to read value for workflow option '$optionName' as a Boolean"
        )
        false
    }
  }

  def getReferenceInputsToMountedPathsOpt(
    createParameters: CreateBatchJobParameters
  ): Option[Map[GcpBatchInput, String]] =
    if (useReferenceDisks) {
      batchAttributes.referenceFileToDiskImageMappingOpt
        .map(getReferenceInputsToMountedPathMappings(_, createParameters.inputOutputParameters.fileInputParameters))
    } else {
      None
    }

  private def generateGcsLocalizationScript(inputs: List[GcpBatchInput],
                                            referenceInputsToMountedPathsOpt: Option[Map[GcpBatchInput, String]]
  )(implicit gcsTransferConfiguration: GcsTransferConfiguration): String = {
    // Generate a mapping of reference inputs to their mounted paths and a section of the localization script to
    // "faux localize" these reference inputs with symlinks to their locations on mounted reference disks.
    import cromwell.backend.google.batch.runnable.RunnableUtils.shellEscaped

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
    outputs: List[GcpBatchOutput]
  )(implicit gcsTransferConfiguration: GcsTransferConfiguration): String = {
    val bundleFunction = (gcsDelocalizationTransferBundle(gcsTransferConfiguration) _).tupled
    generateGcsTransferScript(outputs, bundleFunction) |> bracketTransfersWithMessages("Delocalization")
  }

  private def generateGcsTransferScript[T <: BatchParameter](items: List[T],
                                                             bundleFunction: ((String, NonEmptyList[T])) => String
  ): String = {
    val gcsItems = items collect { case i if i.cloudPath.isInstanceOf[GcsPath] => i }
    groupParametersByGcsBucket(gcsItems) map bundleFunction mkString "\n"
  }

  def uploadGcsLocalizationScript(createParameters: CreateBatchJobParameters,
                                  cloudPath: Path,
                                  transferLibraryContainerPath: Path,
                                  gcsTransferConfiguration: GcsTransferConfiguration,
                                  referenceInputsToMountedPathsOpt: Option[Map[GcpBatchInput, String]]
  ): Future[Unit] = {
    val content = generateGcsLocalizationScript(createParameters.inputOutputParameters.fileInputParameters,
                                                referenceInputsToMountedPathsOpt
    )(gcsTransferConfiguration)
    asyncIo.writeAsync(cloudPath,
                       s"source '$transferLibraryContainerPath'\n\n" + content,
                       Seq(CloudStorageOptions.withMimeType("text/plain"))
    )
  }

  def uploadGcsDelocalizationScript(createParameters: CreateBatchJobParameters,
                                    cloudPath: Path,
                                    transferLibraryContainerPath: Path,
                                    gcsTransferConfiguration: GcsTransferConfiguration
  ): Future[Unit] = {
    val content = generateGcsDelocalizationScript(createParameters.inputOutputParameters.fileOutputParameters)(
      gcsTransferConfiguration
    )
    asyncIo.writeAsync(cloudPath,
                       s"source '$transferLibraryContainerPath'\n\n" + content,
                       Seq(CloudStorageOptions.withMimeType("text/plain"))
    )
  }

  // TAG DISK
  private def createBatchParameters(inputOutputParameters: InputOutputParameters,
                                    customLabels: Seq[GcpLabel]
  ): CreateBatchJobParameters =
    standardParams.backendInitializationDataOption match {
      case Some(data: GcpBackendInitializationData) =>
        val dockerKeyAndToken: Option[CreateBatchDockerKeyAndToken] = for {
          key <- data.privateDockerEncryptionKeyName
          token <- data.privateDockerEncryptedToken
        } yield CreateBatchDockerKeyAndToken(key, token)

        val inputFilePaths = inputOutputParameters.jobInputParameters.map(_.cloudPath.pathAsString).toSet

        val referenceDisksToMount =
          batchAttributes.referenceFileToDiskImageMappingOpt.map(getReferenceDisksToMount(_, inputFilePaths))

        val dockerhubCredentials: (String, String) =
          new String(Base64.getDecoder.decode(batchAttributes.dockerhubToken), "UTF-8").split(":", 2) match {
            case Array(username, password) => (username, password)
            case _ => ("", "")
          }

        val workflowOptions = workflowDescriptor.workflowOptions

        val monitoringImage = new MonitoringImage(
          jobDescriptor = jobDescriptor,
          workflowOptions = workflowOptions,
          workflowPaths = workflowPaths,
          commandDirectory = commandDirectory,
          workingDisk = workingDisk,
          localMonitoringImageScriptPath = localMonitoringImageScriptPath
        )

        val checkpointingConfiguration =
          new CheckpointingConfiguration(
            jobDescriptor = jobDescriptor,
            workflowPaths = workflowPaths,
            commandDirectory = commandDirectory,
            batchConfiguration.batchAttributes.checkpointingInterval
          )

        val enableSshAccess = workflowOptions.getBoolean(WorkflowOptionKeys.EnableSSHAccess).toOption.contains(true)

        // if the `memory_retry_multiplier` is not present in the workflow options there is no need to check whether or
        // not the `stderr` file contained memory retry error keys
        val retryWithMoreMemoryKeys: Option[List[String]] = memoryRetryFactor.flatMap(_ => memoryRetryErrorKeys)

        CreateBatchJobParameters(
          jobDescriptor = jobDescriptor,
          runtimeAttributes = runtimeAttributes,
          dockerImage = jobDockerImage,
          cloudWorkflowRoot = workflowPaths.workflowRoot,
          cloudCallRoot = callRootPath,
          commandScriptContainerPath = cmdInput.containerPath,
          logGcsPath = gcpBatchLogPath,
          inputOutputParameters = inputOutputParameters,
          projectId = googleProject(jobDescriptor.workflowDescriptor),
          computeServiceAccount = computeServiceAccount(jobDescriptor.workflowDescriptor),
          googleLabels = backendLabels ++ customLabels,
          preemptible = preemptible,
          batchTimeout = batchConfiguration.batchTimeout,
          jobShell = batchConfiguration.jobShell,
          privateDockerKeyAndEncryptedToken = dockerKeyAndToken,
          womOutputRuntimeExtractor = jobDescriptor.workflowDescriptor.outputRuntimeExtractor,
          disks = runtimeAttributes.disks,
          virtualPrivateCloudConfiguration = batchAttributes.virtualPrivateCloudConfiguration,
          retryWithMoreMemoryKeys = retryWithMoreMemoryKeys,
          fuseEnabled = fuseEnabled(jobDescriptor.workflowDescriptor),
          referenceDisksForLocalizationOpt = referenceDisksToMount,
          monitoringImage = monitoringImage,
          checkpointingConfiguration,
          enableSshAccess = enableSshAccess,
          vpcNetworkAndSubnetworkProjectLabels = data.vpcNetworkAndSubnetworkProjectLabels,
          dockerhubCredentials = dockerhubCredentials
        )
      case Some(other) =>
        throw new RuntimeException(s"Unexpected initialization data: $other")
      case None =>
        throw new RuntimeException("No batch backend initialization data found?")
    }

  protected def relativePathAndAttachedDisk(path: String,
                                            disks: Seq[GcpBatchAttachedDisk]
  ): (Path, GcpBatchAttachedDisk) = {
    val absolutePath = DefaultPathBuilder.get(path) match {
      case p if !p.isAbsolute => GcpBatchWorkingDisk.MountPoint.resolve(p)
      case p => p
    }

    disks.find(d => absolutePath.startsWith(d.mountPoint)) match {
      case Some(disk) => (disk.mountPoint.relativize(absolutePath), disk)
      case None =>
        throw new Exception(
          s"Absolute path $path doesn't appear to be under any mount points: ${disks.map(_.toString).mkString(", ")}"
        )
    }
  }

  protected def makeSafeReferenceName(referenceName: String): String =
    if (referenceName.length <= 127) referenceName else referenceName.md5Sum

  // De-localize the glob directory as a GcpBatchDirectoryOutput instead of using * pattern match
  protected def generateGlobFileOutputs(womFile: WomGlobFile): List[GcpBatchOutput] = {
    val globName = GlobFunctions.globName(womFile.value)
    val globDirectory = globName + "/"
    val globListFile = globName + ".list"
    val gcsGlobDirectoryDestinationPath = callRootPath.resolve(globDirectory)
    val gcsGlobListFileDestinationPath = callRootPath.resolve(globListFile)

    val (_, globDirectoryDisk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)

    // We need both the glob directory and the glob list:
    List(
      // The glob directory:
      GcpBatchDirectoryOutput(
        makeSafeReferenceName(globDirectory),
        gcsGlobDirectoryDestinationPath,
        DefaultPathBuilder.get(globDirectory),
        globDirectoryDisk,
        optional = false,
        secondary = false
      ),
      // The glob list file:
      GcpBatchFileOutput(
        makeSafeReferenceName(globListFile),
        gcsGlobListFileDestinationPath,
        DefaultPathBuilder.get(globListFile),
        globDirectoryDisk,
        optional = false,
        secondary = false
      )
    )
  }

  lazy val batchMonitoringParamName: String = GcpBatchJobPaths.BatchMonitoringKey
  lazy val localMonitoringLogPath: Path = DefaultPathBuilder.get(gcpBatchCallPaths.batchMonitoringLogFilename)
  lazy val localMonitoringScriptPath: Path = DefaultPathBuilder.get(gcpBatchCallPaths.batchMonitoringScriptFilename)

  lazy val monitoringScript: Option[GcpBatchFileInput] =
    gcpBatchCallPaths.workflowPaths.monitoringScriptPath map { path =>
      GcpBatchFileInput(s"$batchMonitoringParamName-in", path, localMonitoringScriptPath, workingDisk)
    }

  private val DockerMonitoringLogPath: Path =
    GcpBatchWorkingDisk.MountPoint.resolve(gcpBatchCallPaths.batchMonitoringLogFilename)
  private val DockerMonitoringScriptPath: Path =
    GcpBatchWorkingDisk.MountPoint.resolve(gcpBatchCallPaths.batchMonitoringScriptFilename)

  // noinspection ActorMutableStateInspection
  private var hasDockerCredentials: Boolean = false

  override def scriptPreamble: ErrorOr[ScriptPreambleData] =
    if (monitoringOutput.isDefined)
      ScriptPreambleData(s"""|touch $DockerMonitoringLogPath
                             |chmod u+x $DockerMonitoringScriptPath
                             |$DockerMonitoringScriptPath > $DockerMonitoringLogPath &""".stripMargin).valid
    else ScriptPreambleData("").valid

  private[actors] def generateInputs(jobDescriptor: BackendJobDescriptor): Set[GcpBatchInput] = {
    // We need to tell Batch about files that were created as part of command instantiation (these need to be defined
    // as inputs that will be localized down to the VM). Make up 'names' for these files that are just the short
    // md5's of their paths.
    val writeFunctionFiles = instantiatedCommand.createdFiles map { f => f.file.value.md5SumShort -> List(f) } toMap

    val writeFunctionInputs = writeFunctionFiles flatMap { case (name, files) =>
      gcpBatchInputsFromWomFiles(name, files.map(_.file), files.map(localizationPath), jobDescriptor)
    }

    val callInputInputs = callInputFiles flatMap { case (name, files) =>
      gcpBatchInputsFromWomFiles(name, files, files.map(relativeLocalizationPath), jobDescriptor)
    }

    (writeFunctionInputs ++ callInputInputs).toSet
  }

  // Simply create a GcpBatchDirectoryOutput instead of globbing
  protected def generateUnlistedDirectoryOutputs(unlistedDirectory: WomUnlistedDirectory,
                                                 fileEvaluation: FileEvaluation
  ): List[GcpBatchOutput] = {
    val destination = callRootPath.resolve(unlistedDirectory.value.stripPrefix("/"))
    val (relpath, disk) = relativePathAndAttachedDisk(unlistedDirectory.value, runtimeAttributes.disks)
    val directoryOutput = GcpBatchDirectoryOutput(makeSafeReferenceName(unlistedDirectory.value),
                                                  destination,
                                                  relpath,
                                                  disk,
                                                  fileEvaluation.optional,
                                                  fileEvaluation.secondary
    )
    List(directoryOutput)
  }

  private def maybeListedDirectoryToBatchParameters(inputName: String,
                                                    womMaybeListedDirectory: WomMaybeListedDirectory,
                                                    localPath: String
  ) = womMaybeListedDirectory match {
    // If there is a path, simply localize as a directory
    case WomMaybeListedDirectory(Some(path), _, _, _) =>
      List(GcpBatchDirectoryInput(inputName, getPath(path).get, DefaultPathBuilder.get(localPath), workingDisk))

    // If there is a listing, recurse and call gcpBatchInputsFromWomFiles on all the listed files
    case WomMaybeListedDirectory(_, Some(listing), _, _) if listing.nonEmpty =>
      listing.flatMap {
        case womFile: WomFile if isAdHocFile(womFile) =>
          gcpBatchInputsFromWomFiles(makeSafeReferenceName(womFile.valueString),
                                     List(womFile),
                                     List(fileName(womFile)),
                                     jobDescriptor
          )
        case womFile: WomFile =>
          gcpBatchInputsFromWomFiles(makeSafeReferenceName(womFile.valueString),
                                     List(womFile),
                                     List(relativeLocalizationPath(womFile)),
                                     jobDescriptor
          )
      }
    case _ => List.empty
  }

  def generateSingleFileOutputs(womFile: WomSingleFile, fileEvaluation: FileEvaluation): List[GcpBatchFileOutput] = {
    val (relpath, disk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)
    // If the file is on a custom mount point, resolve it so that the full mount path will show up in the cloud path
    // For the default one (cromwell_root), the expectation is that it does not appear
    val mountedPath =
      if (!disk.mountPoint.isSamePathAs(GcpBatchWorkingDisk.Default.mountPoint)) disk.mountPoint.resolve(relpath)
      else relpath
    // Normalize the local path (to get rid of ".." and "."). Also strip any potential leading / so that it gets appended to the call root
    val normalizedPath = mountedPath.normalize().pathAsString.stripPrefix("/")
    val destination = callRootPath.resolve(normalizedPath)
    val batchFileOutput = GcpBatchFileOutput(makeSafeReferenceName(womFile.value),
                                             destination,
                                             relpath,
                                             disk,
                                             fileEvaluation.optional,
                                             fileEvaluation.secondary
    )
    List(batchFileOutput)
  }

  private[actors] def generateOutputs(jobDescriptor: BackendJobDescriptor): Set[GcpBatchOutput] = {
    def evaluateFiles(output: OutputDefinition): List[FileEvaluation] =
      Try(
        output.expression.evaluateFiles(jobDescriptor.localInputs, NoIoFunctionSet, output.womType).map(_.toList)
      ).getOrElse(List.empty[FileEvaluation].validNel)
        .getOrElse(List.empty)

    def relativeFileEvaluation(evaluation: FileEvaluation): FileEvaluation =
      evaluation.copy(file = relativeLocalizationPath(evaluation.file))

    val womFileOutputs = jobDescriptor.taskCall.callable.outputs.flatMap(evaluateFiles) map relativeFileEvaluation

    val outputs: Seq[GcpBatchOutput] = womFileOutputs.distinct flatMap { fileEvaluation =>
      fileEvaluation.file.flattenFiles flatMap {
        case unlistedDirectory: WomUnlistedDirectory =>
          generateUnlistedDirectoryOutputs(unlistedDirectory, fileEvaluation)
        case singleFile: WomSingleFile => generateSingleFileOutputs(singleFile, fileEvaluation)
        case globFile: WomGlobFile => generateGlobFileOutputs(globFile) // Assumes optional = false for globs.
      }
    }

    val additionalGlobOutput =
      jobDescriptor.taskCall.callable.additionalGlob.toList.flatMap(generateGlobFileOutputs).toSet

    outputs.toSet ++ additionalGlobOutput
  }

  protected def uploadGcsTransferLibrary(createBatchParameters: CreateBatchJobParameters,
                                         cloudPath: Path,
                                         gcsTransferConfiguration: GcsTransferConfiguration
  ): Future[Unit] =
    asyncIo.writeAsync(cloudPath, gcsTransferLibrary, Seq(CloudStorageOptions.withMimeType("text/plain")))

  lazy val monitoringOutput: Option[GcpBatchFileOutput] = monitoringScript map { _ =>
    GcpBatchFileOutput(
      s"$batchMonitoringParamName-out",
      gcpBatchCallPaths.batchMonitoringLogPath,
      localMonitoringLogPath,
      workingDisk,
      optional = false,
      secondary = false,
      contentType = plainTextContentType
    )
  }

  override lazy val commandDirectory: Path = GcpBatchWorkingDisk.MountPoint

  // Primary entry point for cromwell to run GCP Batch job
  override def executeAsync(): Future[ExecutionHandle] = {

    // Want to force runtimeAttributes to evaluate so we can fail quickly now if we need to:
    def evaluateRuntimeAttributes = Future.fromTry(Try(runtimeAttributes))

    def generateInputOutputParameters: Future[InputOutputParameters] = Future.fromTry(Try {
      val rcFileOutput = GcpBatchFileOutput(
        returnCodeFilename,
        returnCodeGcsPath,
        DefaultPathBuilder.get(returnCodeFilename),
        workingDisk,
        optional = false,
        secondary = false,
        contentType = plainTextContentType
      )

      val memoryRetryRCFileOutput = GcpBatchFileOutput(
        memoryRetryRCFilename,
        memoryRetryRCGcsPath,
        DefaultPathBuilder.get(memoryRetryRCFilename),
        workingDisk,
        optional = true,
        secondary = false,
        contentType = plainTextContentType
      )

      case class StandardStream(name: String, f: StandardPaths => Path) {
        val filename: String = f(gcpBatchCallPaths.standardPaths).name
      }

      val standardStreams = List(
        StandardStream("stdout", _.output),
        StandardStream("stderr", _.error)
      ) map { s =>
        GcpBatchFileOutput(
          s.name,
          returnCodeGcsPath.sibling(s.filename),
          DefaultPathBuilder.get(s.filename),
          workingDisk,
          optional = false,
          secondary = false,
          uploadPeriod = batchAttributes.logFlushPeriod,
          contentType = plainTextContentType
        )
      }

      InputOutputParameters(
        DetritusInputParameters(
          executionScriptInputParameter = cmdInput,
          monitoringScriptInputParameter = monitoringScript
        ),
        generateInputs(jobDescriptor).toList,
        standardStreams ++ generateOutputs(jobDescriptor).toList,
        DetritusOutputParameters(
          monitoringScriptOutputParameter = monitoringOutput,
          rcFileOutputParameter = rcFileOutput,
          memoryRetryRCFileOutputParameter = memoryRetryRCFileOutput
        ),
        List.empty
      )

    })

    val gcpBatchParameters = CreateGcpBatchParameters(
      jobDescriptor = jobDescriptor,
      runtimeAttributes = runtimeAttributes,
      batchAttributes = batchAttributes,
      projectId = batchAttributes.project,
      region = batchAttributes.location
    )

    val runBatchResponse = for {
      _ <- evaluateRuntimeAttributes
      _ <- uploadScriptFile()
      customLabels <- Future.fromTry(GcpLabel.fromWorkflowOptions(workflowDescriptor.workflowOptions))
      batchParameters <- generateInputOutputParameters
      createParameters = createBatchParameters(batchParameters, customLabels)
      drsLocalizationManifestCloudPath = jobPaths.callExecutionRoot / GcpBatchJobPaths.DrsLocalizationManifestName
      _ <- uploadDrsLocalizationManifest(createParameters, drsLocalizationManifestCloudPath)
      gcsTransferConfiguration = initializationData.gcpBatchConfiguration.batchAttributes.gcsTransferConfiguration
      gcsTransferLibraryCloudPath = jobPaths.callExecutionRoot / GcpBatchJobPaths.GcsTransferLibraryName
      transferLibraryContainerPath = createParameters.commandScriptContainerPath.sibling(GcsTransferLibraryName)
      _ <- uploadGcsTransferLibrary(createParameters, gcsTransferLibraryCloudPath, gcsTransferConfiguration)
      gcsLocalizationScriptCloudPath = jobPaths.callExecutionRoot / GcpBatchJobPaths.GcsLocalizationScriptName
      referenceInputsToMountedPathsOpt = getReferenceInputsToMountedPathsOpt(createParameters)
      _ <- uploadGcsLocalizationScript(createParameters,
                                       gcsLocalizationScriptCloudPath,
                                       transferLibraryContainerPath,
                                       gcsTransferConfiguration,
                                       referenceInputsToMountedPathsOpt
      )
      gcsDelocalizationScriptCloudPath = jobPaths.callExecutionRoot / GcpBatchJobPaths.GcsDelocalizationScriptName
      _ <- uploadGcsDelocalizationScript(createParameters,
                                         gcsDelocalizationScriptCloudPath,
                                         transferLibraryContainerPath,
                                         gcsTransferConfiguration
      )
      _ = this.hasDockerCredentials = createParameters.privateDockerKeyAndEncryptedToken.isDefined
      jobName = "job-" + java.util.UUID.randomUUID.toString
      request = GcpBatchRequest(workflowId, createParameters, jobName = jobName, gcpBatchParameters)
      response <- runBatchJob(request = request,
                              backendSingletonActor = backendSingletonActor,
                              requestFactory = initializationData.requestFactory,
                              jobLogger = jobLogger
      )
      _ = sendGoogleLabelsToMetadata(customLabels)
      _ = sendIncrementMetricsForReferenceFiles(referenceInputsToMountedPathsOpt.map(_.keySet))

    } yield response

    runBatchResponse
      .map { runId =>
        PendingExecutionHandle(jobDescriptor = jobDescriptor,
                               pendingJob = runId,
                               runInfo = Option(Run(runId)),
                               previousState = None
        )
      }
      .recover { case JobAbortedException => AbortedExecutionHandle }
  }

  val futureKvJobKey: KvJobKey =
    KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt + 1)

  override def reconnectAsync(jobId: StandardAsyncJob): Future[ExecutionHandle] = {
    log.info("reconnect async runs") // in for debugging remove later
    val handle = PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState](
      jobDescriptor,
      jobId,
      Option(Run(jobId)),
      previousState = None
    )
    Future.successful(handle)
  }

  override lazy val pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(5.second, 5.minutes, 1.1)

  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff =
    SimpleExponentialBackoff(initialInterval = 5.seconds, maxInterval = 20.seconds, multiplier = 1.1)

  protected def sendIncrementMetricsForReferenceFiles(referenceInputFilesOpt: Option[Set[GcpBatchInput]]): Unit =
    referenceInputFilesOpt match {
      case Some(referenceInputFiles) =>
        referenceInputFiles.foreach { referenceInputFile =>
          increment(NonEmptyList.of("referencefiles", referenceInputFile.relativeHostPath.pathAsString))
        }
      case _ =>
      // do nothing - reference disks feature is either not configured in Cromwell or disabled in workflow options
    }

  protected def sendIncrementMetricsForDockerImageCache(dockerImageCacheDiskOpt: Option[String],
                                                        dockerImageAsSpecifiedByUser: String,
                                                        isDockerImageCacheUsageRequested: Boolean
  ): Unit =
    (isDockerImageCacheUsageRequested, dockerImageCacheDiskOpt) match {
      case (true, None) =>
        increment(NonEmptyList("docker", List("image", "cache", "image_not_in_cache", dockerImageAsSpecifiedByUser)))
      case (true, Some(_)) =>
        increment(NonEmptyList("docker", List("image", "cache", "used_image_from_cache", dockerImageAsSpecifiedByUser)))
      case (false, Some(_)) =>
        increment(NonEmptyList("docker", List("image", "cache", "cached_image_not_used", dockerImageAsSpecifiedByUser)))
      case _ => // docker image cache not requested and image is not in cache anyway - do nothing
    }

  override def pollStatusAsync(handle: GcpBatchPendingExecutionHandle): Future[RunStatus] = {
    // yes, we use the whole jobName as the id
    val jobNameStr = handle.pendingJob.jobId

    for {
      _ <- Future.unit // trick to get into a future context
      _ = log.info(s"started polling for $jobNameStr")
      jobName = JobName.parse(jobNameStr)
      status <- pollStatus(workflowId, jobName, backendSingletonActor, initializationData.requestFactory)
    } yield status // RunStatus.fromJobStatus(job.getStatus.getState)
  }

  override def isTerminal(runStatus: RunStatus): Boolean =
    runStatus match {
      case _: RunStatus.TerminalRunStatus =>
        log.info(s"isTerminal match terminal run status with $runStatus")
        true
      case other =>
        log.info(f"isTerminal match _ running with status $other")
        false
    }

  override def isDone(runStatus: RunStatus): Boolean =
    runStatus match {
      case _: RunStatus.Success => true
      case _: RunStatus.UnsuccessfulRunStatus => false
      case _ =>
        throw new RuntimeException(
          s"Cromwell programmer blunder: isSuccess was called on an incomplete RunStatus ($runStatus)."
        )
    }
//  override def isDone(runStatus: RunStatus): Boolean =
//    runStatus match {
//      case _: RunStatus.Succeeded =>
//        log.info("GCP batch job succeeded matched isDone")
//        true
//      case _: RunStatus.UnsuccessfulRunStatus =>
//        log.info("GCP batch job unsuccessful matched isDone")
//        false
//      case _ =>
//        log.info(s"did not match isDone: $runStatus")
//        throw new RuntimeException(
//          s"Cromwell programmer blunder: isDone was called on an incomplete RunStatus ($runStatus)."
//        )
//    }

  override def getTerminalEvents(runStatus: RunStatus): Seq[ExecutionEvent] =
    runStatus match {
      case successStatus: RunStatus.Success => successStatus.eventList
      case unknown =>
        throw new RuntimeException(s"handleExecutionSuccess not called with RunStatus.Success. Instead got $unknown")
    }
//  override def getTerminalEvents(runStatus: RunStatus): Seq[ExecutionEvent] =
//    runStatus match {
//      case t: RunStatus.TerminalRunStatus =>
//        log.warning(s"Tried to get terminal events on a terminal status without events: $runStatus")
//        t.eventList
//
//      case unknown =>
//        throw new RuntimeException(s"handleExecutionSuccess not called with RunStatus.Success. Instead got $unknown")
//    }

  // TODO: Find methods that are overriden by papi but we are not overriding them
  override def retryEvaluateOutputs(exception: Exception): Boolean =
    exception match {
      case aggregated: CromwellAggregatedException =>
        aggregated.throwables.collectFirst { case s: SocketTimeoutException => s }.isDefined
      case _ => false
    }

  private lazy val standardPaths = jobPaths.standardPaths

  override def handleExecutionFailure(runStatus: RunStatus, returnCode: Option[Int]): Future[ExecutionHandle] = {

    def generateBetterErrorMsg(runStatus: RunStatus.UnsuccessfulRunStatus, errorMsg: String): String =
      if (
        runStatus.errorCode.getCode.value() == BatchFailedPreConditionErrorCode
        && errorMsg.contains("Execution failed")
        && (errorMsg.contains("Localization") || errorMsg.contains("Delocalization"))
      ) {
        s"Please check the log file for more details: $gcpBatchLogPath."
      }
      // If error code 10, add some extra messaging to the server logging
      else if (runStatus.errorCode.getCode.value() == BatchMysteriouslyCrashedErrorCode) {
        jobLogger.info(s"Job Failed with Error Code 10 for a machine where Preemptible is set to $preemptible")
        errorMsg
      } else errorMsg

    // Inner function: Handles a 'Failed' runStatus (or Preempted if preemptible was false)
    def handleFailedRunStatus(runStatus: RunStatus.UnsuccessfulRunStatus, returnCode: Option[Int]): ExecutionHandle = {

      lazy val prettyError = runStatus.prettyPrintedError

      def isDockerPullFailure: Boolean = prettyError.contains("not found: does not exist or no pull access")

      (runStatus.errorCode, runStatus.jesCode) match {
        case (GrpcStatus.NOT_FOUND, Some(JesFailedToDelocalize)) =>
          FailedNonRetryableExecutionHandle(
            FailedToDelocalizeFailure(prettyError, jobTag, Option(standardPaths.error)),
            kvPairsToSave = None
          )
        case (GrpcStatus.ABORTED, Some(JesUnexpectedTermination)) =>
          handleUnexpectedTermination(runStatus.errorCode, prettyError, returnCode)
        case _ if isDockerPullFailure =>
          val unable = s"Unable to pull Docker image '$jobDockerImage' "
          val details =
            if (hasDockerCredentials)
              "but Docker credentials are present; is this Docker account authorized to pull the image? "
            else
              "and there are effectively no Docker credentials present (one or more of token, authorization, or Google KMS key may be missing). " +
                "Please check your private Docker configuration and/or the pull access for this image. "
          val message = unable + details + prettyError
          FailedNonRetryableExecutionHandle(
            StandardException(runStatus.errorCode, message, jobTag, returnCode, standardPaths.error),
            returnCode,
            None
          )
        case _ =>
          val finalPrettyPrintedError = generateBetterErrorMsg(runStatus, prettyError)
          FailedNonRetryableExecutionHandle(
            StandardException(runStatus.errorCode, finalPrettyPrintedError, jobTag, returnCode, standardPaths.error),
            returnCode,
            None
          )
      }
    }

    Future.fromTry {
      Try {
        runStatus match {
          case preemptedStatus: RunStatus.Preempted if preemptible => handlePreemption(preemptedStatus, returnCode)
          case _: RunStatus.Cancelled => AbortedExecutionHandle
          case failedStatus: RunStatus.UnsuccessfulRunStatus => handleFailedRunStatus(failedStatus, returnCode)
          case unknown =>
            throw new RuntimeException(
              s"handleExecutionFailure not called with RunStatus.Failed or RunStatus.Preempted. Instead got $unknown"
            )
        }
      }
    }
  }

  private def nextAttemptPreemptedAndUnexpectedRetryCountsToKvPairs(p: Int, ur: Int): Seq[KvPair] =
    Seq(
      KvPair(ScopedKey(workflowId, futureKvJobKey, GcpBatchBackendLifecycleActorFactory.unexpectedRetryCountKey),
             ur.toString
      ),
      KvPair(ScopedKey(workflowId, futureKvJobKey, GcpBatchBackendLifecycleActorFactory.preemptionCountKey), p.toString)
    )

  private def handleUnexpectedTermination(errorCode: GrpcStatus,
                                          errorMessage: String,
                                          jobReturnCode: Option[Int]
  ): ExecutionHandle = {
    val msg = s"Retrying. $errorMessage"
    previousRetryReasons match {
      case Valid(PreviousRetryReasons(p, ur)) =>
        val thisUnexpectedRetry = ur + 1
        if (thisUnexpectedRetry <= maxUnexpectedRetries) {
          val preemptionAndUnexpectedRetryCountsKvPairs =
            nextAttemptPreemptedAndUnexpectedRetryCountsToKvPairs(p, thisUnexpectedRetry)
          // Increment unexpected retry count and preemption count stays the same
          FailedRetryableExecutionHandle(
            StandardException(errorCode, msg, jobTag, jobReturnCode, standardPaths.error),
            jobReturnCode,
            kvPairsToSave = Option(preemptionAndUnexpectedRetryCountsKvPairs)
          )
        } else {
          FailedNonRetryableExecutionHandle(
            StandardException(errorCode, errorMessage, jobTag, jobReturnCode, standardPaths.error),
            jobReturnCode,
            None
          )
        }
      case Invalid(_) =>
        FailedNonRetryableExecutionHandle(
          StandardException(errorCode, errorMessage, jobTag, jobReturnCode, standardPaths.error),
          jobReturnCode,
          None
        )
    }
  }

  private def handlePreemption(runStatus: RunStatus.Preempted, jobReturnCode: Option[Int]): ExecutionHandle = {
    import common.numeric.IntegerUtil._

    val errorCode: GrpcStatus = runStatus.errorCode
    val prettyPrintedError: String = runStatus.prettyPrintedError
    previousRetryReasons match {
      case Valid(PreviousRetryReasons(p, ur)) =>
        val thisPreemption = p + 1
        val taskName = s"${workflowDescriptor.id}:${call.localName}"
        val baseMsg = s"Task $taskName was preempted for the ${thisPreemption.toOrdinal} time."

        val preemptionAndUnexpectedRetryCountsKvPairs =
          nextAttemptPreemptedAndUnexpectedRetryCountsToKvPairs(thisPreemption, ur)
        if (thisPreemption < maxPreemption) {
          // Increment preemption count and unexpectedRetryCount stays the same
          val msg =
            s"$baseMsg The call will be restarted with another preemptible VM (max preemptible attempts number is " +
              s"$maxPreemption). Error code $errorCode.$prettyPrintedError"
          FailedRetryableExecutionHandle(
            StandardException(errorCode, msg, jobTag, jobReturnCode, standardPaths.error),
            jobReturnCode,
            kvPairsToSave = Option(preemptionAndUnexpectedRetryCountsKvPairs)
          )
        } else {
          val msg = s"$baseMsg The maximum number of preemptible attempts ($maxPreemption) has been reached. The " +
            s"call will be restarted with a non-preemptible VM. Error code $errorCode.$prettyPrintedError)"
          FailedRetryableExecutionHandle(
            StandardException(errorCode, msg, jobTag, jobReturnCode, standardPaths.error),
            jobReturnCode,
            kvPairsToSave = Option(preemptionAndUnexpectedRetryCountsKvPairs)
          )
        }
      case Invalid(_) =>
        FailedNonRetryableExecutionHandle(
          StandardException(errorCode, prettyPrintedError, jobTag, jobReturnCode, standardPaths.error),
          jobReturnCode,
          None
        )
    }
  }

  override lazy val startMetadataKeyValues: Map[String, Any] =
    super[GcpBatchJobCachingActorHelper].startMetadataKeyValues

  // TODO: review sending machine and Instance type
  override def getTerminalMetadata(runStatus: RunStatus): Map[String, Any] =
    runStatus match {
      case _: TerminalRunStatus => Map()
      case unknown => throw new RuntimeException(s"Attempt to get terminal metadata from non terminal status: $unknown")
    }

  override def mapOutputWomFile(womFile: WomFile): WomFile =
    womFileToGcsPath(generateOutputs(jobDescriptor))(womFile)

  override def globParentDirectory(womGlobFile: WomGlobFile): Path = {
    val (_, disk) = relativePathAndAttachedDisk(womGlobFile.value, runtimeAttributes.disks)
    disk.mountPoint
  }

  protected def googleProject(descriptor: BackendWorkflowDescriptor): String =
    descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleProject, batchAttributes.project)

  protected def computeServiceAccount(descriptor: BackendWorkflowDescriptor): String =
    descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleComputeServiceAccount,
                                         batchAttributes.computeServiceAccount
    )

  protected def fuseEnabled(descriptor: BackendWorkflowDescriptor): Boolean =
    descriptor.workflowOptions.getBoolean(WorkflowOptionKeys.EnableFuse).toOption.getOrElse(batchAttributes.enableFuse)

  protected def useDockerImageCache(descriptor: BackendWorkflowDescriptor): Boolean =
    descriptor.workflowOptions.getBoolean(WorkflowOptionKeys.UseDockerImageCache).getOrElse(false)

  override def cloudResolveWomFile(womFile: WomFile): WomFile =
    womFile.mapFile { value =>
      getPath(value) match {
        case Success(drsPath: DrsPath) => DrsResolver.getSimpleGsUri(drsPath).unsafeRunSync().getOrElse(value)
        case Success(path) => path.pathAsString
        case _ => value
      }
    }

  override def mapCommandLineWomFile(womFile: WomFile): WomFile =
    womFile.mapFile { value =>
      (getPath(value), asAdHocFile(womFile)) match {
        case (Success(gcsPath: GcsPath), Some(adHocFile)) =>
          // Ad hoc files will be placed directly at the root ("/cromwell_root/ad_hoc_file.txt") unlike other input files
          // for which the full path is being propagated ("/cromwell_root/path/to/input_file.txt")
          workingDisk.mountPoint.resolve(adHocFile.alternativeName.getOrElse(gcsPath.name)).pathAsString
        case (Success(path @ (_: GcsPath | _: HttpPath)), _) =>
          workingDisk.mountPoint.resolve(path.pathWithoutScheme).pathAsString
        case (Success(drsPath: DrsPath), _) =>
          val filePath = DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
          workingDisk.mountPoint.resolve(filePath).pathAsString
        case (Success(sraPath: SraPath), _) =>
          workingDisk.mountPoint.resolve(s"sra-${sraPath.accession}/${sraPath.pathWithoutScheme}").pathAsString
        case _ => value
      }
    }

  override def mapCommandLineJobInputWomFile(womFile: WomFile): WomFile =
    womFile.mapFile(value =>
      getPath(value) match {
        case Success(gcsPath: GcsPath) => workingDisk.mountPoint.resolve(gcsPath.pathWithoutScheme).pathAsString
        case Success(drsPath: DrsPath) =>
          val filePath = DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
          workingDisk.mountPoint.resolve(filePath).pathAsString
        case _ => value
      }
    )

  def womFileToGcsPath(batchOutputs: Set[GcpBatchOutput])(womFile: WomFile): WomFile =
    womFile mapFile { path =>
      batchOutputs collectFirst {
        case batchOutput if batchOutput.name == makeSafeReferenceName(path) =>
          val pathAsString = batchOutput.cloudPath.pathAsString
          if (batchOutput.isFileParameter && !batchOutput.cloudPath.exists) {
            // This is not an error if the path represents a `File?` optional output (the Batch delocalization script
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

  // No need for Cromwell-performed localization in the Batch backend, ad hoc values are localized directly from GCS to the VM by Batch.
  override lazy val localizeAdHocValues: List[AdHocValue] => ErrorOr[List[StandardAdHocValue]] =
    _.map(Coproduct[StandardAdHocValue](_)).validNel
}
