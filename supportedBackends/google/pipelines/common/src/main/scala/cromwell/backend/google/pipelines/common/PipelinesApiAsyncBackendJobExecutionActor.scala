package cromwell.backend.google.pipelines.common

import java.net.SocketTimeoutException
import _root_.io.grpc.Status
import akka.actor.ActorRef
import akka.http.scaladsl.model.{ContentType, ContentTypes}
import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import cats.syntax.validated._
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.cloud.storage.contrib.nio.CloudStorageOptions
import common.util.StringUtil._
import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.backend._
import cromwell.backend.async.{
  AbortedExecutionHandle,
  ExecutionHandle,
  FailedNonRetryableExecutionHandle,
  FailedRetryableExecutionHandle,
  PendingExecutionHandle
}
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common.PipelinesApiJobPaths.GcsTransferLibraryName
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory._
import cromwell.backend.google.pipelines.common.api.RunStatus.TerminalRunStatus
import cromwell.backend.google.pipelines.common.api._
import cromwell.backend.google.pipelines.common.api.clients.PipelinesApiRunCreationClient.JobAbortedException
import cromwell.backend.google.pipelines.common.api.clients.{
  PipelinesApiAbortClient,
  PipelinesApiRunCreationClient,
  PipelinesApiStatusRequestClient
}
import cromwell.backend.google.pipelines.common.authentication.PipelinesApiDockerCredentials
import cromwell.backend.google.pipelines.common.errors.FailedToDelocalizeFailure
import cromwell.backend.google.pipelines.common.io._
import cromwell.backend.google.pipelines.common.monitoring.{CheckpointingConfiguration, MonitoringImage}
import cromwell.backend.io.DirectoryFunctions
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
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import cromwell.filesystems.http.HttpPath
import cromwell.filesystems.sra.SraPath
import cromwell.google.pipelines.common.PreviousRetryReasons
import cromwell.services.keyvalue.KeyValueServiceActor._
import cromwell.services.metadata.CallMetadataKeys
import mouse.all._
import shapeless.Coproduct
import wom.callable.Callable.OutputDefinition
import wom.callable.MetaValueElement.{MetaValueElementBoolean, MetaValueElementObject}
import wom.callable.{AdHocValue, RuntimeEnvironment}
import wom.core.FullyQualifiedName
import wom.expression.{FileEvaluation, NoIoFunctionSet}
import wom.types.{WomArrayType, WomSingleFileType}
import wom.values._

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object PipelinesApiAsyncBackendJobExecutionActor {
  val JesOperationIdKey = "__jes_operation_id"

  type JesPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, Run, RunStatus]

  val maxUnexpectedRetries = 2

  val JesFailedToDelocalize = 5
  val JesUnexpectedTermination = 13
  val JesPreemption = 14

  val PapiFailedPreConditionErrorCode = 9
  val PapiMysteriouslyCrashedErrorCode = 10

  // If the JES code is 2 (UNKNOWN), this sub-string indicates preemption:
  val FailedToStartDueToPreemptionSubstring = "failed to start due to preemption"
  val FailedV2Style = "The assigned worker has failed to complete the operation"

  val plainTextContentType: Option[ContentType.WithCharset] = Option(ContentTypes.`text/plain(UTF-8)`)

  def StandardException(errorCode: Status,
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

    new Exception(s"Task $jobTag failed. $returnCodeMessage PAPI error code ${errorCode.getCode.value}. $message")
  }
}

class PipelinesApiAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
    extends BackendJobLifecycleActor
    with StandardAsyncExecutionActor
    with PipelinesApiJobCachingActorHelper
    with PipelinesApiStatusRequestClient
    with PipelinesApiRunCreationClient
    with PipelinesApiAbortClient
    with PipelinesApiReferenceFilesMappingOperations
    with PipelinesApiDockerCacheMappingOperations
    with PapiInstrumentation {

  override lazy val ioCommandBuilder: IoCommandBuilder = GcsBatchCommandBuilder

  import PipelinesApiAsyncBackendJobExecutionActor._

  override lazy val workflowId: WorkflowId = jobDescriptor.workflowDescriptor.id
  override lazy val requestFactory: PipelinesApiRequestFactory = initializationData.genomicsRequestFactory

  val jesBackendSingletonActor: ActorRef =
    standardParams.backendSingletonActorOption.getOrElse(
      throw new RuntimeException("JES Backend actor cannot exist without the JES backend singleton actor")
    )

  override type StandardAsyncRunInfo = Run

  override type StandardAsyncRunState = RunStatus

  def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean = thiz == that

  override val papiApiActor: ActorRef = jesBackendSingletonActor

  override lazy val pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(
    initialInterval = 30 seconds,
    maxInterval = jesAttributes.maxPollingInterval seconds,
    multiplier = 1.1
  )

  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff =
    SimpleExponentialBackoff(initialInterval = 3 seconds, maxInterval = 20 seconds, multiplier = 1.1)

  override lazy val runtimeEnvironment: RuntimeEnvironment =
    RuntimeEnvironmentBuilder(jobDescriptor.runtimeAttributes,
                              PipelinesApiWorkingDisk.MountPoint,
                              PipelinesApiWorkingDisk.MountPoint
    )(standardParams.minimumRuntimeSettings)

  protected lazy val cmdInput: PipelinesApiFileInput =
    PipelinesApiFileInput(PipelinesApiJobPaths.JesExecParamName,
                          pipelinesApiCallPaths.script,
                          DefaultPathBuilder.get(pipelinesApiCallPaths.scriptFilename),
                          workingDisk
    )

  protected lazy val dockerConfiguration: Option[PipelinesApiDockerCredentials] =
    pipelinesConfiguration.dockerCredentials

  protected val previousRetryReasons: ErrorOr[PreviousRetryReasons] =
    PreviousRetryReasons.tryApply(jobDescriptor.prefetchedKvStoreEntries, jobDescriptor.key.attempt)

  protected lazy val jobDockerImage: String =
    jobDescriptor.maybeCallCachingEligible.dockerHash.getOrElse(runtimeAttributes.dockerImage)

  override lazy val dockerImageUsed: Option[String] = Option(jobDockerImage)

  override lazy val preemptible: Boolean = previousRetryReasons match {
    case Valid(PreviousRetryReasons(p, _)) => p < maxPreemption
    case _ => false
  }

  override def tryAbort(job: StandardAsyncJob): Unit = abortJob(job)

  override def requestsAbortAndDiesImmediately: Boolean = false

  /*_*/ // Silence an errant IntelliJ warning
  override def receive: Receive =
    pollingActorClientReceive orElse runCreationClientReceive orElse abortActorClientReceive orElse kvClientReceive orElse super.receive
  /*_*/ // https://stackoverflow.com/questions/36679973/controlling-false-intellij-code-editor-error-in-scala-plugin

  /**
    * Takes two arrays of remote and local WOM File paths and generates the necessary `PipelinesApiInput`s.
    */
  protected def pipelinesApiInputsFromWomFiles(jesNamePrefix: String,
                                               remotePathArray: Seq[WomFile],
                                               localPathArray: Seq[WomFile],
                                               jobDescriptor: BackendJobDescriptor
  ): Iterable[PipelinesApiInput] =
    (remotePathArray zip localPathArray zipWithIndex) flatMap { case ((remotePath, localPath), index) =>
      Seq(
        PipelinesApiFileInput(s"$jesNamePrefix-$index",
                              getPath(remotePath.valueString).get,
                              DefaultPathBuilder.get(localPath.valueString),
                              workingDisk
        )
      )
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

  override protected def fileName(file: WomFile): WomFile =
    file.mapFile(value =>
      getPath(value) match {
        case Success(drsPath: DrsPath) =>
          DefaultPathBuilder.get(DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()).name
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

  protected def callInputFiles: Map[FullyQualifiedName, Seq[WomFile]] =
    jobDescriptor.fullyQualifiedInputs map { case (key, womFile) =>
      val arrays: Seq[WomArray] = womFile collectAsSeq {
        case womFile: WomFile if !inputsToNotLocalize.contains(womFile) =>
          val files: List[WomSingleFile] = DirectoryFunctions
            .listWomSingleFiles(womFile, pipelinesApiCallPaths.workflowPaths)
            .toTry(s"Error getting single files for $womFile")
            .get
          WomArray(WomArrayType(WomSingleFileType), files)
      }

      key -> arrays.flatMap(_.value).collect { case womFile: WomFile =>
        womFile
      }
    }

  private[pipelines] def generateInputs(jobDescriptor: BackendJobDescriptor): Set[PipelinesApiInput] = {
    // We need to tell PAPI about files that were created as part of command instantiation (these need to be defined
    // as inputs that will be localized down to the VM). Make up 'names' for these files that are just the short
    // md5's of their paths.
    val writeFunctionFiles = instantiatedCommand.createdFiles map { f => f.file.value.md5SumShort -> List(f) } toMap

    val writeFunctionInputs = writeFunctionFiles flatMap { case (name, files) =>
      pipelinesApiInputsFromWomFiles(name, files.map(_.file), files.map(localizationPath), jobDescriptor)
    }

    val callInputInputs = callInputFiles flatMap { case (name, files) =>
      pipelinesApiInputsFromWomFiles(name, files, files.map(relativeLocalizationPath), jobDescriptor)
    }

    (writeFunctionInputs ++ callInputInputs).toSet
  }

  /**
    * Given a path (relative or absolute), returns a (Path, JesAttachedDisk) tuple where the Path is
    * relative to the AttachedDisk's mount point
    *
    * @throws Exception if the `path` does not live in one of the supplied `disks`
    */
  protected def relativePathAndAttachedDisk(path: String,
                                            disks: Seq[PipelinesApiAttachedDisk]
  ): (Path, PipelinesApiAttachedDisk) = {
    val absolutePath = DefaultPathBuilder.get(path) match {
      case p if !p.isAbsolute => PipelinesApiWorkingDisk.MountPoint.resolve(p)
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

  /**
    * If the desired reference name is too long, we don't want to break JES or risk collisions by arbitrary truncation. So,
    * just use a hash. We only do this when needed to give better traceability in the normal case.
    */
  protected def makeSafeReferenceName(referenceName: String): String =
    if (referenceName.length <= 127) referenceName else referenceName.md5Sum

  protected[pipelines] def generateOutputs(jobDescriptor: BackendJobDescriptor): Set[PipelinesApiOutput] = {
    import cats.syntax.validated._
    def evaluateFiles(output: OutputDefinition): List[FileEvaluation] =
      Try(
        output.expression.evaluateFiles(jobDescriptor.localInputs, NoIoFunctionSet, output.womType).map(_.toList)
      ).getOrElse(List.empty[FileEvaluation].validNel)
        .getOrElse(List.empty)

    def relativeFileEvaluation(evaluation: FileEvaluation): FileEvaluation =
      evaluation.copy(file = relativeLocalizationPath(evaluation.file))

    val womFileOutputs = jobDescriptor.taskCall.callable.outputs.flatMap(evaluateFiles) map relativeFileEvaluation

    val outputs: Seq[PipelinesApiOutput] = womFileOutputs.distinct flatMap { fileEvaluation =>
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

  protected def generateUnlistedDirectoryOutputs(womFile: WomUnlistedDirectory,
                                                 fileEvaluation: FileEvaluation
  ): List[PipelinesApiOutput] = {
    val directoryPath = womFile.value.ensureSlashed
    val directoryListFile = womFile.value.ensureUnslashed + ".list"
    val gcsDirDestinationPath = callRootPath.resolve(directoryPath)
    val gcsListDestinationPath = callRootPath.resolve(directoryListFile)

    val (_, directoryDisk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)

    // We need both the collection directory and the collection list:
    List(
      // The collection directory:
      PipelinesApiFileOutput(
        makeSafeReferenceName(directoryListFile),
        gcsListDestinationPath,
        DefaultPathBuilder.get(directoryListFile),
        directoryDisk,
        fileEvaluation.optional,
        fileEvaluation.secondary
      ),
      // The collection list file:
      PipelinesApiFileOutput(
        makeSafeReferenceName(directoryPath),
        gcsDirDestinationPath,
        DefaultPathBuilder.get(directoryPath + "*"),
        directoryDisk,
        fileEvaluation.optional,
        fileEvaluation.secondary
      )
    )
  }

  protected def generateSingleFileOutputs(womFile: WomSingleFile,
                                          fileEvaluation: FileEvaluation
  ): List[PipelinesApiFileOutput] = {
    val destination = callRootPath.resolve(womFile.value.stripPrefix("/"))
    val (relpath, disk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)
    val jesFileOutput = PipelinesApiFileOutput(makeSafeReferenceName(womFile.value),
                                               destination,
                                               relpath,
                                               disk,
                                               fileEvaluation.optional,
                                               fileEvaluation.secondary
    )
    List(jesFileOutput)
  }

  protected def generateGlobFileOutputs(womFile: WomGlobFile): List[PipelinesApiOutput] = {
    val globName = GlobFunctions.globName(womFile.value)
    val globDirectory = globName + "/"
    val globListFile = globName + ".list"
    val gcsGlobDirectoryDestinationPath = callRootPath.resolve(globDirectory)
    val gcsGlobListFileDestinationPath = callRootPath.resolve(globListFile)

    val (_, globDirectoryDisk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)

    // We need both the glob directory and the glob list:
    List(
      // The glob directory:
      PipelinesApiFileOutput(
        makeSafeReferenceName(globDirectory),
        gcsGlobDirectoryDestinationPath,
        DefaultPathBuilder.get(globDirectory + "*"),
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

  lazy val jesMonitoringParamName: String = PipelinesApiJobPaths.JesMonitoringKey
  lazy val localMonitoringLogPath: Path = DefaultPathBuilder.get(pipelinesApiCallPaths.jesMonitoringLogFilename)
  lazy val localMonitoringScriptPath: Path = DefaultPathBuilder.get(pipelinesApiCallPaths.jesMonitoringScriptFilename)
  lazy val localMonitoringImageScriptPath: Path =
    DefaultPathBuilder.get(pipelinesApiCallPaths.jesMonitoringImageScriptFilename)

  lazy val monitoringScript: Option[PipelinesApiFileInput] =
    pipelinesApiCallPaths.workflowPaths.monitoringScriptPath map { path =>
      PipelinesApiFileInput(s"$jesMonitoringParamName-in", path, localMonitoringScriptPath, workingDisk)
    }

  lazy val monitoringOutput: Option[PipelinesApiFileOutput] = monitoringScript map { _ =>
    PipelinesApiFileOutput(
      s"$jesMonitoringParamName-out",
      pipelinesApiCallPaths.jesMonitoringLogPath,
      localMonitoringLogPath,
      workingDisk,
      optional = false,
      secondary = false,
      contentType = plainTextContentType
    )
  }

  override lazy val commandDirectory: Path = PipelinesApiWorkingDisk.MountPoint

  private val DockerMonitoringLogPath: Path =
    PipelinesApiWorkingDisk.MountPoint.resolve(pipelinesApiCallPaths.jesMonitoringLogFilename)
  private val DockerMonitoringScriptPath: Path =
    PipelinesApiWorkingDisk.MountPoint.resolve(pipelinesApiCallPaths.jesMonitoringScriptFilename)
  // noinspection ActorMutableStateInspection
  private var hasDockerCredentials: Boolean = false

  private lazy val isDockerImageCacheUsageRequested =
    runtimeAttributes.useDockerImageCache.getOrElse(useDockerImageCache(jobDescriptor.workflowDescriptor))

  override def scriptPreamble: ErrorOr[ScriptPreambleData] =
    if (monitoringOutput.isDefined)
      ScriptPreambleData(s"""|touch $DockerMonitoringLogPath
                             |chmod u+x $DockerMonitoringScriptPath
                             |$DockerMonitoringScriptPath > $DockerMonitoringLogPath &""".stripMargin).valid
    else ScriptPreambleData("").valid

  override def globParentDirectory(womGlobFile: WomGlobFile): Path = {
    val (_, disk) = relativePathAndAttachedDisk(womGlobFile.value, runtimeAttributes.disks)
    disk.mountPoint
  }

  protected def googleProject(descriptor: BackendWorkflowDescriptor): String =
    descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleProject, jesAttributes.project)

  protected def computeServiceAccount(descriptor: BackendWorkflowDescriptor): String =
    descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleComputeServiceAccount,
                                         jesAttributes.computeServiceAccount
    )

  protected def fuseEnabled(descriptor: BackendWorkflowDescriptor): Boolean =
    descriptor.workflowOptions.getBoolean(WorkflowOptionKeys.EnableFuse).toOption.getOrElse(jesAttributes.enableFuse)

  protected def googleLegacyMachineSelection(descriptor: BackendWorkflowDescriptor): Boolean =
    descriptor.workflowOptions.getBoolean(WorkflowOptionKeys.GoogleLegacyMachineSelection).getOrElse(false)

  protected def useDockerImageCache(descriptor: BackendWorkflowDescriptor): Boolean =
    descriptor.workflowOptions.getBoolean(WorkflowOptionKeys.UseDockerImageCache).getOrElse(false)

  override def isTerminal(runStatus: RunStatus): Boolean =
    runStatus match {
      case _: TerminalRunStatus => true
      case _ => false
    }

  private def createPipelineParameters(inputOutputParameters: InputOutputParameters,
                                       customLabels: Seq[GoogleLabel]
  ): CreatePipelineParameters =
    standardParams.backendInitializationDataOption match {
      case Some(data: PipelinesApiBackendInitializationData) =>
        val dockerKeyAndToken: Option[CreatePipelineDockerKeyAndToken] = for {
          key <- data.privateDockerEncryptionKeyName
          token <- data.privateDockerEncryptedToken
        } yield CreatePipelineDockerKeyAndToken(key, token)

        val inputFilePaths = inputOutputParameters.jobInputParameters.map(_.cloudPath.pathAsString).toSet
        val referenceDisksToMount =
          if (useReferenceDisks)
            jesAttributes.referenceFileToDiskImageMappingOpt.map(getReferenceDisksToMount(_, inputFilePaths))
          else None

        val workflowOptions = workflowDescriptor.workflowOptions

        val monitoringImage =
          new MonitoringImage(
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
            pipelinesConfiguration.papiAttributes.checkpointingInterval
          )

        val enableSshAccess = workflowOptions.getBoolean(WorkflowOptionKeys.EnableSSHAccess).toOption.contains(true)

        val dockerImageCacheDiskOpt =
          getDockerCacheDiskImageForAJob(
            dockerImageToCacheDiskImageMappingOpt = jesAttributes.dockerImageToCacheDiskImageMappingOpt,
            dockerImageAsSpecifiedByUser = runtimeAttributes.dockerImage,
            dockerImageWithDigest = jobDockerImage,
            jobLogger = jobLogger
          )

        // if the `memory_retry_multiplier` is not present in the workflow options there is no need to check whether or
        // not the `stderr` file contained memory retry error keys
        val retryWithMoreMemoryKeys: Option[List[String]] = memoryRetryFactor.flatMap(_ => memoryRetryErrorKeys)

        CreatePipelineParameters(
          jobDescriptor = jobDescriptor,
          runtimeAttributes = runtimeAttributes,
          dockerImage = jobDockerImage,
          cloudWorkflowRoot = workflowPaths.workflowRoot,
          cloudCallRoot = callRootPath,
          commandScriptContainerPath = cmdInput.containerPath,
          logGcsPath = jesLogPath,
          inputOutputParameters = inputOutputParameters,
          projectId = googleProject(jobDescriptor.workflowDescriptor),
          computeServiceAccount = computeServiceAccount(jobDescriptor.workflowDescriptor),
          googleLabels = backendLabels ++ customLabels,
          preemptible = preemptible,
          pipelineTimeout = pipelinesConfiguration.pipelineTimeout,
          jobShell = pipelinesConfiguration.jobShell,
          privateDockerKeyAndEncryptedToken = dockerKeyAndToken,
          womOutputRuntimeExtractor = jobDescriptor.workflowDescriptor.outputRuntimeExtractor,
          disks = runtimeAttributes.disks,
          virtualPrivateCloudConfiguration = jesAttributes.virtualPrivateCloudConfiguration,
          retryWithMoreMemoryKeys = retryWithMoreMemoryKeys,
          fuseEnabled = fuseEnabled(jobDescriptor.workflowDescriptor),
          referenceDisksForLocalizationOpt = referenceDisksToMount,
          monitoringImage = monitoringImage,
          checkpointingConfiguration,
          enableSshAccess = enableSshAccess,
          vpcNetworkAndSubnetworkProjectLabels = data.vpcNetworkAndSubnetworkProjectLabels,
          dockerImageCacheDiskOpt = isDockerImageCacheUsageRequested.option(dockerImageCacheDiskOpt).flatten
        )
      case Some(other) =>
        throw new RuntimeException(s"Unexpected initialization data: $other")
      case None =>
        throw new RuntimeException("No pipelines API backend initialization data found?")
    }

  override def isFatal(throwable: Throwable): Boolean = super.isFatal(throwable) || isFatalJesException(throwable)

  override def isTransient(throwable: Throwable): Boolean = isTransientJesException(throwable)

  override def executeAsync(): Future[ExecutionHandle] = createNewJob()

  val futureKvJobKey: KvJobKey =
    KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt + 1)

  override def recoverAsync(jobId: StandardAsyncJob): Future[ExecutionHandle] = reconnectToExistingJob(jobId)

  override def reconnectAsync(jobId: StandardAsyncJob): Future[ExecutionHandle] = reconnectToExistingJob(jobId)

  override def reconnectToAbortAsync(jobId: StandardAsyncJob): Future[ExecutionHandle] =
    reconnectToExistingJob(jobId, forceAbort = true)

  private def reconnectToExistingJob(jobForResumption: StandardAsyncJob, forceAbort: Boolean = false) = {
    if (forceAbort) tryAbort(jobForResumption)
    Future.successful(
      PendingExecutionHandle(jobDescriptor, jobForResumption, Option(Run(jobForResumption)), previousState = None)
    )
  }

  protected def uploadDrsLocalizationManifest(createPipelineParameters: CreatePipelineParameters,
                                              cloudPath: Path
  ): Future[Unit] = Future.successful(())

  protected def uploadGcsTransferLibrary(createPipelineParameters: CreatePipelineParameters,
                                         cloudPath: Path,
                                         gcsTransferConfiguration: GcsTransferConfiguration
  ): Future[Unit] = Future.successful(())

  protected def uploadGcsLocalizationScript(createPipelineParameters: CreatePipelineParameters,
                                            cloudPath: Path,
                                            transferLibraryContainerPath: Path,
                                            gcsTransferConfiguration: GcsTransferConfiguration,
                                            referenceInputsToMountedPathsOpt: Option[Map[PipelinesApiInput, String]]
  ): Future[Unit] = Future.successful(())

  protected def uploadGcsDelocalizationScript(createPipelineParameters: CreatePipelineParameters,
                                              cloudPath: Path,
                                              transferLibraryContainerPath: Path,
                                              gcsTransferConfiguration: GcsTransferConfiguration
  ): Future[Unit] = Future.successful(())

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

  private def createNewJob(): Future[ExecutionHandle] = {
    // Want to force runtimeAttributes to evaluate so we can fail quickly now if we need to:
    def evaluateRuntimeAttributes = Future.fromTry(Try(runtimeAttributes))

    def generateInputOutputParameters: Future[InputOutputParameters] = Future.fromTry(Try {
      val rcFileOutput = PipelinesApiFileOutput(
        returnCodeFilename,
        returnCodeGcsPath,
        DefaultPathBuilder.get(returnCodeFilename),
        workingDisk,
        optional = false,
        secondary = false,
        contentType = plainTextContentType
      )

      val memoryRetryRCFileOutput = PipelinesApiFileOutput(
        memoryRetryRCFilename,
        memoryRetryRCGcsPath,
        DefaultPathBuilder.get(memoryRetryRCFilename),
        workingDisk,
        optional = true,
        secondary = false,
        contentType = plainTextContentType
      )

      case class StandardStream(name: String, f: StandardPaths => Path) {
        val filename: String = f(pipelinesApiCallPaths.standardPaths).name
      }

      val standardStreams = List(
        StandardStream("stdout", _.output),
        StandardStream("stderr", _.error)
      ) map { s =>
        PipelinesApiFileOutput(
          s.name,
          returnCodeGcsPath.sibling(s.filename),
          DefaultPathBuilder.get(s.filename),
          workingDisk,
          optional = false,
          secondary = false,
          uploadPeriod = jesAttributes.logFlushPeriod,
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

    def uploadScriptFile = commandScriptContents.fold(
      errors => Future.failed(new RuntimeException(errors.toList.mkString(", "))),
      asyncIo.writeAsync(jobPaths.script, _, Seq(CloudStorageOptions.withMimeType("text/plain")))
    )

    def sendGoogleLabelsToMetadata(customLabels: Seq[GoogleLabel]): Unit = {
      lazy val backendLabelEvents: Map[String, String] =
        ((backendLabels ++ customLabels) map { l => s"${CallMetadataKeys.BackendLabels}:${l.key}" -> l.value }).toMap
      tellMetadata(backendLabelEvents)
    }

    def getReferenceInputsToMountedPathsOpt(
      createPipelinesParameters: CreatePipelineParameters
    ): Option[Map[PipelinesApiInput, String]] =
      if (useReferenceDisks) {
        jesAttributes.referenceFileToDiskImageMappingOpt
          .map(
            getReferenceInputsToMountedPathMappings(_,
                                                    createPipelinesParameters.inputOutputParameters.fileInputParameters
            )
          )
      } else {
        None
      }

    val runPipelineResponse = for {
      _ <- evaluateRuntimeAttributes
      _ <- uploadScriptFile
      customLabels <- Future.fromTry(GoogleLabels.fromWorkflowOptions(workflowDescriptor.workflowOptions))
      jesParameters <- generateInputOutputParameters
      createParameters = createPipelineParameters(jesParameters, customLabels)
      drsLocalizationManifestCloudPath = jobPaths.callExecutionRoot / PipelinesApiJobPaths.DrsLocalizationManifestName
      _ <- uploadDrsLocalizationManifest(createParameters, drsLocalizationManifestCloudPath)
      gcsTransferConfiguration = initializationData.papiConfiguration.papiAttributes.gcsTransferConfiguration
      gcsTransferLibraryCloudPath = jobPaths.callExecutionRoot / PipelinesApiJobPaths.GcsTransferLibraryName
      transferLibraryContainerPath = createParameters.commandScriptContainerPath.sibling(GcsTransferLibraryName)
      _ <- uploadGcsTransferLibrary(createParameters, gcsTransferLibraryCloudPath, gcsTransferConfiguration)
      gcsLocalizationScriptCloudPath = jobPaths.callExecutionRoot / PipelinesApiJobPaths.GcsLocalizationScriptName
      referenceInputsToMountedPathsOpt = getReferenceInputsToMountedPathsOpt(createParameters)
      _ <- uploadGcsLocalizationScript(createParameters,
                                       gcsLocalizationScriptCloudPath,
                                       transferLibraryContainerPath,
                                       gcsTransferConfiguration,
                                       referenceInputsToMountedPathsOpt
      )
      gcsDelocalizationScriptCloudPath = jobPaths.callExecutionRoot / PipelinesApiJobPaths.GcsDelocalizationScriptName
      _ <- uploadGcsDelocalizationScript(createParameters,
                                         gcsDelocalizationScriptCloudPath,
                                         transferLibraryContainerPath,
                                         gcsTransferConfiguration
      )
      _ = this.hasDockerCredentials = createParameters.privateDockerKeyAndEncryptedToken.isDefined
      runId <- runPipeline(workflowId, createParameters, jobLogger)
      _ = sendGoogleLabelsToMetadata(customLabels)
      _ = sendIncrementMetricsForReferenceFiles(referenceInputsToMountedPathsOpt.map(_.keySet))
      _ = sendIncrementMetricsForDockerImageCache(
        dockerImageCacheDiskOpt = createParameters.dockerImageCacheDiskOpt,
        dockerImageAsSpecifiedByUser = runtimeAttributes.dockerImage,
        isDockerImageCacheUsageRequested = isDockerImageCacheUsageRequested
      )
    } yield runId

    runPipelineResponse map { runId =>
      PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None)
    } recover { case JobAbortedException =>
      AbortedExecutionHandle
    }
  }

  protected def sendIncrementMetricsForReferenceFiles(referenceInputFilesOpt: Option[Set[PipelinesApiInput]]): Unit =
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

  override def pollStatusAsync(handle: JesPendingExecutionHandle): Future[RunStatus] =
    super[PipelinesApiStatusRequestClient].pollStatus(workflowId, handle.pendingJob)

  override def customPollStatusFailure: PartialFunction[(ExecutionHandle, Exception), ExecutionHandle] = {
    case (_: JesPendingExecutionHandle @unchecked, JobAbortedException) =>
      AbortedExecutionHandle
    case (oldHandle: JesPendingExecutionHandle @unchecked, e: GoogleJsonResponseException) if e.getStatusCode == 404 =>
      jobLogger.error(s"JES Job ID ${oldHandle.runInfo.get.job} has not been found, failing call")
      FailedNonRetryableExecutionHandle(e, kvPairsToSave = None)
  }

  override lazy val startMetadataKeyValues: Map[String, Any] =
    super[PipelinesApiJobCachingActorHelper].startMetadataKeyValues

  override def getTerminalMetadata(runStatus: RunStatus): Map[String, Any] =
    runStatus match {
      case terminalRunStatus: TerminalRunStatus =>
        Map(
          PipelinesApiMetadataKeys.MachineType -> terminalRunStatus.machineType.getOrElse("unknown"),
          PipelinesApiMetadataKeys.InstanceName -> terminalRunStatus.instanceName.getOrElse("unknown"),
          PipelinesApiMetadataKeys.Zone -> terminalRunStatus.zone.getOrElse("unknown")
        )
      case unknown => throw new RuntimeException(s"Attempt to get terminal metadata from non terminal status: $unknown")
    }

  override def mapOutputWomFile(womFile: WomFile): WomFile =
    womFileToGcsPath(generateOutputs(jobDescriptor))(womFile)

  protected[pipelines] def womFileToGcsPath(jesOutputs: Set[PipelinesApiOutput])(womFile: WomFile): WomFile =
    womFile mapFile { path =>
      jesOutputs collectFirst {
        case jesOutput if jesOutput.name == makeSafeReferenceName(path) => jesOutput.cloudPath.pathAsString
      } getOrElse path
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

  override def getTerminalEvents(runStatus: RunStatus): Seq[ExecutionEvent] =
    runStatus match {
      case successStatus: RunStatus.Success => successStatus.eventList
      case unknown =>
        throw new RuntimeException(s"handleExecutionSuccess not called with RunStatus.Success. Instead got $unknown")
    }

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
        runStatus.errorCode.getCode.value == PapiFailedPreConditionErrorCode
        && errorMsg.contains("Execution failed")
        && (errorMsg.contains("Localization") || errorMsg.contains("Delocalization"))
      ) {
        s"Please check the log file for more details: $jesLogPath."
      }
      // If error code 10, add some extra messaging to the server logging
      else if (runStatus.errorCode.getCode.value == PapiMysteriouslyCrashedErrorCode) {
        jobLogger.info(s"Job Failed with Error Code 10 for a machine where Preemptible is set to $preemptible")
        errorMsg
      } else errorMsg

    // Inner function: Handles a 'Failed' runStatus (or Preempted if preemptible was false)
    def handleFailedRunStatus(runStatus: RunStatus.UnsuccessfulRunStatus, returnCode: Option[Int]): ExecutionHandle = {

      lazy val prettyError = runStatus.prettyPrintedError

      def isDockerPullFailure: Boolean = prettyError.contains("not found: does not exist or no pull access")

      (runStatus.errorCode, runStatus.jesCode) match {
        case (Status.NOT_FOUND, Some(JesFailedToDelocalize)) =>
          FailedNonRetryableExecutionHandle(
            FailedToDelocalizeFailure(prettyError, jobTag, Option(standardPaths.error)),
            kvPairsToSave = None
          )
        case (Status.ABORTED, Some(JesUnexpectedTermination)) =>
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
      KvPair(ScopedKey(workflowId, futureKvJobKey, PipelinesApiBackendLifecycleActorFactory.unexpectedRetryCountKey),
             ur.toString
      ),
      KvPair(ScopedKey(workflowId, futureKvJobKey, PipelinesApiBackendLifecycleActorFactory.preemptionCountKey),
             p.toString
      )
    )

  private def handleUnexpectedTermination(errorCode: Status,
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

    val errorCode: Status = runStatus.errorCode
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

  // No need for Cromwell-performed localization in the PAPI backend, ad hoc values are localized directly from GCS to the VM by PAPI.
  override lazy val localizeAdHocValues: List[AdHocValue] => ErrorOr[List[StandardAdHocValue]] =
    _.map(Coproduct[StandardAdHocValue](_)).validNel
}
