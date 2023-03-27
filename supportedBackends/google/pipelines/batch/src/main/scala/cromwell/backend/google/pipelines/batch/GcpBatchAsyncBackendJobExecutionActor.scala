package cromwell.backend.google.pipelines.batch

import akka.util.Timeout
import akka.http.scaladsl.model.{ContentType, ContentTypes}
//import cats.syntax.validated._
import com.google.api.gax.rpc.NotFoundException
import common.util.StringUtil._
//import common.validation.ErrorOr._
import common.validation.Validation._
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.backend._
import cromwell.core._
import cromwell.backend.google.pipelines.batch.io._
import cromwell.backend.google.pipelines.batch.GcpBatchJobPaths.GcsTransferLibraryName
import cats.data.Validated.Valid
//import common.validation.Validation._
import cromwell.backend.google.pipelines.batch.GcpBatchRequestFactory._
import cromwell.backend.google.pipelines.common.monitoring.CheckpointingConfiguration
//import cromwell.backend.google.pipelines.common.monitoring.{CheckpointingConfiguration, MonitoringImage}
import java.util.concurrent.ExecutionException
import scala.concurrent.Await
import cromwell.core.{ExecutionEvent, WorkflowId}
import cromwell.backend.async.PendingExecutionHandle
import cromwell.backend.async.ExecutionHandle
import akka.actor.ActorRef
import akka.pattern.AskSupport
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.backend.google.pipelines.batch.GcpBatchConfigurationAttributes.GcsTransferConfiguration
import scala.concurrent.Future
import scala.concurrent.duration._
import GcpBatchBackendSingletonActor._
import cats.data.NonEmptyList
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.google.pipelines.batch.RunStatus.{Running, Succeeded, TerminalRunStatus}
import cromwell.backend.google.pipelines.common.WorkflowOptionKeys
import cromwell.core.io.IoCommandBuilder
import cromwell.core.path.DefaultPathBuilder
import cromwell.filesystems.drs.{DrsPath, DrsResolver}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wom.values.WomFile
import wom.types.{WomArrayType, WomSingleFileType}
import wom.values._
import wom.core.FullyQualifiedName
import cromwell.backend.io.DirectoryFunctions
import wom.callable.MetaValueElement.{MetaValueElementBoolean, MetaValueElementObject}
import cromwell.core.path.Path

import scala.util.{Success, Try}
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.http.HttpPath
import cromwell.filesystems.sra.SraPath

import wom.format.MemorySize
import wom.expression.{FileEvaluation, NoIoFunctionSet}
import wom.callable.Callable.OutputDefinition
import wdl4s.parser.MemoryUnit

import scala.language.postfixOps

import cromwell.backend.google.pipelines.common.GoogleLabels

object GcpBatchAsyncBackendJobExecutionActor {
  val JesOperationIdKey = "__jes_operation_id"

  type GcpBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, Run, RunStatus]

  val plainTextContentType: Option[ContentType.WithCharset] = Option(ContentTypes.`text/plain(UTF-8)`)

}

class GcpBatchAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends BackendJobLifecycleActor
    with StandardAsyncExecutionActor
    with AskSupport
    with GcpBatchJobCachingActorHelper
    with GcpBatchStatusRequestClient
    with CromwellInstrumentation {

  import GcpBatchAsyncBackendJobExecutionActor._
  override lazy val ioCommandBuilder: IoCommandBuilder = GcsBatchCommandBuilder

  lazy val gcpBatchCommand: String = jobDescriptor.taskCall.callable.commandTemplateString(Map.empty)
  lazy val workflowId: WorkflowId = jobDescriptor.workflowDescriptor.id
  //lazy val gcpBootDiskSizeGb = runtimeAttributes.bootDiskSize
  //println(f"$gcpBootDiskSizeGb LOOOOOOOOOOOK")


  /** The type of the run info when a job is started. */
  override type StandardAsyncRunInfo = Run

  /** The type of the run status returned during each poll. */
  override type StandardAsyncRunState = RunStatus

  // temporary until GCP Batch client can generate random job IDs
  private val jobTemp = "job-" + java.util.UUID.randomUUID.toString

  //override def receive: Receive = pollingActorClientReceive orElse runCreationClientReceive orElse super.receive
  //override def receive: Receive = pollingActorClientReceive orElse super.receive

  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */
  def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean = thiz == that

  protected lazy val cmdInput: GcpBatchFileInput =
    GcpBatchFileInput(GcpBatchJobPaths.BatchExecParamName, gcpBatchCallPaths.script, DefaultPathBuilder.get(gcpBatchCallPaths.scriptFilename), workingDisk)

  private lazy val jobDockerImage = jobDescriptor.maybeCallCachingEligible.dockerHash
                                                 .getOrElse(runtimeAttributes.dockerImage)

  override def dockerImageUsed: Option[String] = Option(jobDockerImage)

  // Need to add previousRetryReasons and preemptible in order to get preemptible to work in the tests
  protected val previousRetryReasons: ErrorOr[PreviousRetryReasons] = PreviousRetryReasons.tryApply(jobDescriptor.prefetchedKvStoreEntries, jobDescriptor.key.attempt)


   lazy val preemptible: Boolean = previousRetryReasons match {
    case Valid(PreviousRetryReasons(p, _)) => p < maxPreemption
    case _ => false
  }
  //private var hasDockerCredentials: Boolean = false

  //type GcpBatchPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, Run, StandardAsyncRunState]

  val backendSingletonActor: ActorRef = standardParams.backendSingletonActorOption
                                                      .getOrElse(throw new RuntimeException("GCP Batch actor cannot exist without its backend singleton 2"))


  /**
   * Takes two arrays of remote and local WOM File paths and generates the necessary `PipelinesApiInput`s.
   */
  protected def gcpBatchInputsFromWomFiles(jesNamePrefix: String,
                                           remotePathArray: Seq[WomFile],
                                           localPathArray: Seq[WomFile],
                                           jobDescriptor: BackendJobDescriptor): Iterable[GcpBatchInput] = {
    (remotePathArray zip localPathArray zipWithIndex) flatMap {
      case ((remotePath, localPath), index) =>
        Seq(GcpBatchFileInput(s"$jesNamePrefix-$index", getPath(remotePath.valueString).get, DefaultPathBuilder.get(localPath.valueString), workingDisk))
    }
  }


  /**
    * Turns WomFiles into relative paths.  These paths are relative to the working disk.
    *
    * relativeLocalizationPath("foo/bar.txt") -> "foo/bar.txt"
    * relativeLocalizationPath("gs://some/bucket/foo.txt") -> "some/bucket/foo.txt"
    */
  override protected def relativeLocalizationPath(file: WomFile): WomFile = {
    file.mapFile(value =>
      getPath(value) match {
        case Success(drsPath: DrsPath) => DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
        case Success(path) => path.pathWithoutScheme
        case _ => value
      }
    )
  }

  lazy val localMonitoringImageScriptPath: Path =
    DefaultPathBuilder.get(gcpBatchCallPaths.batchMonitoringImageScriptFilename)


  override protected def fileName(file: WomFile): WomFile = {
    file.mapFile(value =>
      getPath(value) match {
        case Success(drsPath: DrsPath) => DefaultPathBuilder
          .get(DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()).name
        case Success(path) => path.name
        case _ => value
      }
    )
  }

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

  protected def callInputFiles: Map[FullyQualifiedName, Seq[WomFile]] = {

    jobDescriptor.fullyQualifiedInputs map {
      case (key, womFile) =>
        val arrays: Seq[WomArray] = womFile collectAsSeq {
          case womFile: WomFile if !inputsToNotLocalize.contains(womFile) =>
            val files: List[WomSingleFile] = DirectoryFunctions
              .listWomSingleFiles(womFile, gcpBatchCallPaths.workflowPaths)
              .toTry(s"Error getting single files for $womFile").get
            WomArray(WomArrayType(WomSingleFileType), files)
        }

        key -> arrays.flatMap(_.value).collect {
          case womFile: WomFile => womFile
        }
    }
  }

  def uploadScriptFile(): Future[Unit] = {
  commandScriptContents
    .fold(
      errors => Future
        .failed(new RuntimeException(errors
          .toList
          .mkString(", "))),
      asyncIo
        .writeAsync(jobPaths
          .script, _, Seq
          .empty)
    )
  }


  private def createPipelineParameters(inputOutputParameters: InputOutputParameters,
                                       //customLabels: Seq[GoogleLabel],
                                      ): CreatePipelineParameters = {
    standardParams.backendInitializationDataOption match {
      case Some(data: GcpBackendInitializationData) =>
        val dockerKeyAndToken: Option[CreatePipelineDockerKeyAndToken] = for {
          key <- data.privateDockerEncryptionKeyName
          token <- data.privateDockerEncryptedToken
        } yield CreatePipelineDockerKeyAndToken(key, token)

        /*
           * Right now this doesn't cost anything, because sizeOption returns the size if it was previously already fetched
           * for some reason (expression evaluation for instance), but otherwise does not retrieve it and returns None.
           * In CWL-land we tend to be aggressive in pre-fetching the size in order to be able to evaluate JS expressions,
           * but less in WDL as we can get it last minute and on demand because size is a WDL function, whereas in CWL
           * we don't inspect the JS to know if size is called and therefore always pre-fetch it.
           *
           * We could decide to call withSize before in which case we would retrieve the size for all files and have
           * a guaranteed more accurate total size, but there might be performance impacts ?
           */
        val inputFileSize = Option(callInputFiles.values.flatMap(_.flatMap(_.sizeOption)).sum)

        // Attempt to adjust the disk size by taking into account the size of input files
        val adjustedSizeDisks = inputFileSize.map(size => MemorySize.apply(size.toDouble, MemoryUnit.Bytes).to(MemoryUnit.GB)) map { inputFileSizeInformation =>
          runtimeAttributes.disks.adjustWorkingDiskWithNewMin(
            inputFileSizeInformation,
            jobLogger.info(s"Adjusted working disk size to ${inputFileSizeInformation.amount} GB to account for input files")
          )
        } getOrElse runtimeAttributes.disks

        //val inputFilePaths = inputOutputParameters.jobInputParameters.map(_.cloudPath.pathAsString).toSet
        /*
        val referenceDisksToMount =
          batchAttributes.referenceFileToDiskImageMappingOpt.map(getReferenceDisksToMount(_, inputFilePaths))

        */

        val workflowOptions = workflowDescriptor.workflowOptions


        /*
        val monitoringImage =
          new MonitoringImage(
            jobDescriptor = jobDescriptor,
            workflowOptions = workflowOptions,
            workflowPaths = workflowPaths,
            commandDirectory = commandDirectory,
            workingDisk = workingDisk,
            localMonitoringImageScriptPath = localMonitoringImageScriptPath,
          )

        */

        val checkpointingConfiguration =
          new CheckpointingConfiguration(
            jobDescriptor = jobDescriptor,
            workflowPaths = workflowPaths,
            commandDirectory = commandDirectory,
            batchConfiguration.batchAttributes.checkpointingInterval
          )

        val enableSshAccess = workflowOptions.getBoolean(WorkflowOptionKeys.EnableSSHAccess).toOption.contains(true)

        /*
        val dockerImageCacheDiskOpt =
          getDockerCacheDiskImageForAJob(
            dockerImageToCacheDiskImageMappingOpt = jesAttributes.dockerImageToCacheDiskImageMappingOpt,
            dockerImageAsSpecifiedByUser = runtimeAttributes.dockerImage,
            dockerImageWithDigest = jobDockerImage,
            jobLogger = jobLogger
          )

        */

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
          logGcsPath = gcpBatchLogPath,
          inputOutputParameters = inputOutputParameters,
          projectId = googleProject(jobDescriptor.workflowDescriptor),
          computeServiceAccount = computeServiceAccount(jobDescriptor.workflowDescriptor),
          //googleLabels = backendLabels ++ customLabels,
          preemptible = preemptible,
          pipelineTimeout = batchConfiguration.pipelineTimeout,
          jobShell = batchConfiguration.jobShell,
          privateDockerKeyAndEncryptedToken = dockerKeyAndToken,
          womOutputRuntimeExtractor = jobDescriptor.workflowDescriptor.outputRuntimeExtractor,
          adjustedSizeDisks = adjustedSizeDisks,
          virtualPrivateCloudConfiguration = batchAttributes.virtualPrivateCloudConfiguration,
          retryWithMoreMemoryKeys = retryWithMoreMemoryKeys,
          fuseEnabled = fuseEnabled(jobDescriptor.workflowDescriptor),
          //referenceDisksForLocalizationOpt = referenceDisksToMount,
          //monitoringImage = monitoringImage,
          checkpointingConfiguration,
          enableSshAccess = enableSshAccess,
          //vpcNetworkAndSubnetworkProjectLabels = data.vpcNetworkAndSubnetworkProjectLabels,
          //dockerImageCacheDiskOpt = isDockerImageCacheUsageRequested.option(dockerImageCacheDiskOpt).flatten
        )
      case Some(other) =>
        throw new RuntimeException(s"Unexpected initialization data: $other")
      case None =>
        throw new RuntimeException("No pipelines API backend initialization data found?")
    }
  }

  protected def relativePathAndAttachedDisk(path: String, disks: Seq[GcpBatchAttachedDisk]): (Path, GcpBatchAttachedDisk) = {
    val absolutePath = DefaultPathBuilder.get(path) match {
      case p if !p.isAbsolute => GcpBatchWorkingDisk.MountPoint.resolve(p)
      case p => p
    }

    disks.find(d => absolutePath.startsWith(d.mountPoint)) match {
      case Some(disk) => (disk.mountPoint.relativize(absolutePath), disk)
      case None =>
        throw new Exception(s"Absolute path $path doesn't appear to be under any mount points: ${disks.map(_.toString).mkString(", ")}")
    }
  }

  protected def makeSafeReferenceName(referenceName: String): String = {
    if (referenceName.length <= 127) referenceName else referenceName.md5Sum
  }

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
      GcpBatchFileOutput(makeSafeReferenceName(globDirectory), gcsGlobDirectoryDestinationPath, DefaultPathBuilder.get(globDirectory + "*"), globDirectoryDisk,
        optional = false, secondary = false),
      // The glob list file:
      GcpBatchFileOutput(makeSafeReferenceName(globListFile), gcsGlobListFileDestinationPath, DefaultPathBuilder.get(globListFile), globDirectoryDisk,
        optional = false, secondary = false)
    )
  }

  lazy val batchMonitoringParamName: String = GcpBatchJobPaths.BatchMonitoringKey
  lazy val localMonitoringLogPath: Path = DefaultPathBuilder.get(gcpBatchCallPaths.batchMonitoringLogFilename)
  lazy val localMonitoringScriptPath: Path = DefaultPathBuilder.get(gcpBatchCallPaths.batchMonitoringScriptFilename)

  lazy val monitoringScript: Option[GcpBatchFileInput] = {
    gcpBatchCallPaths.workflowPaths.monitoringScriptPath map { path =>
      GcpBatchFileInput(s"$batchMonitoringParamName-in", path, localMonitoringScriptPath, workingDisk)
    }
  }

  private def generateInputs(jobDescriptor: BackendJobDescriptor): Set[GcpBatchInput] = {
    // We need to tell PAPI about files that were created as part of command instantiation (these need to be defined
    // as inputs that will be localized down to the VM). Make up 'names' for these files that are just the short
    // md5's of their paths.
    val writeFunctionFiles = instantiatedCommand.createdFiles map { f => f.file.value.md5SumShort -> List(f) } toMap

    val writeFunctionInputs = writeFunctionFiles flatMap {
      case (name, files) => gcpBatchInputsFromWomFiles(name, files.map(_.file), files.map(localizationPath), jobDescriptor)
    }

    val callInputInputs = callInputFiles flatMap {
      case (name, files) => gcpBatchInputsFromWomFiles(name, files, files.map(relativeLocalizationPath), jobDescriptor)
    }

    (writeFunctionInputs ++ callInputInputs).toSet
  }

  protected def generateUnlistedDirectoryOutputs(womFile: WomUnlistedDirectory, fileEvaluation: FileEvaluation): List[GcpBatchOutput] = {
    val directoryPath = womFile.value.ensureSlashed
    val directoryListFile = womFile.value.ensureUnslashed + ".list"
    val gcsDirDestinationPath = callRootPath.resolve(directoryPath)
    val gcsListDestinationPath = callRootPath.resolve(directoryListFile)

    val (_, directoryDisk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)

    // We need both the collection directory and the collection list:
    List(
      // The collection directory:
      GcpBatchFileOutput(
        makeSafeReferenceName(directoryListFile),
        gcsListDestinationPath,
        DefaultPathBuilder.get(directoryListFile),
        directoryDisk,
        fileEvaluation.optional,
        fileEvaluation.secondary
      ),
      // The collection list file:
      GcpBatchFileOutput(
        makeSafeReferenceName(directoryPath),
        gcsDirDestinationPath,
        DefaultPathBuilder.get(directoryPath + "*"),
        directoryDisk,
        fileEvaluation.optional,
        fileEvaluation.secondary
      )
    )
  }

  protected def generateSingleFileOutputs(womFile: WomSingleFile, fileEvaluation: FileEvaluation): List[GcpBatchFileOutput] = {
    val destination = callRootPath.resolve(womFile.value.stripPrefix("/"))
    val (relpath, disk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)
    val jesFileOutput = GcpBatchFileOutput(makeSafeReferenceName(womFile.value), destination, relpath, disk, fileEvaluation.optional, fileEvaluation.secondary)
    List(jesFileOutput)
  }

  protected def generateOutputs(jobDescriptor: BackendJobDescriptor): Set[GcpBatchOutput] = {
    import cats.syntax.validated._
    def evaluateFiles(output: OutputDefinition): List[FileEvaluation] = {
      Try(
        output.expression.evaluateFiles(jobDescriptor.localInputs, NoIoFunctionSet, output.womType).map(_.toList)
      ).getOrElse(List.empty[FileEvaluation].validNel)
        .getOrElse(List.empty)
    }

    def relativeFileEvaluation(evaluation: FileEvaluation): FileEvaluation = {
      evaluation.copy(file = relativeLocalizationPath(evaluation.file))
    }

    val womFileOutputs = jobDescriptor.taskCall.callable.outputs.flatMap(evaluateFiles) map relativeFileEvaluation

    val outputs: Seq[GcpBatchOutput] = womFileOutputs.distinct flatMap { fileEvaluation =>
      fileEvaluation.file.flattenFiles flatMap {
        case unlistedDirectory: WomUnlistedDirectory => generateUnlistedDirectoryOutputs(unlistedDirectory, fileEvaluation)
        case singleFile: WomSingleFile => generateSingleFileOutputs(singleFile, fileEvaluation)
        case globFile: WomGlobFile => generateGlobFileOutputs(globFile) // Assumes optional = false for globs.
      }
    }

    val additionalGlobOutput = jobDescriptor.taskCall.callable.additionalGlob.toList.flatMap(generateGlobFileOutputs).toSet

    outputs.toSet ++ additionalGlobOutput
  }
  protected def uploadGcsTransferLibrary(createPipelineParameters: CreatePipelineParameters, cloudPath: Path, gcsTransferConfiguration: GcsTransferConfiguration): Future[Unit] = Future.successful(())


  //private lazy val standardPaths = jobPaths.standardPaths

  lazy val monitoringOutput: Option[GcpBatchFileOutput] = monitoringScript map { _ =>
    GcpBatchFileOutput(s"$batchMonitoringParamName-out",
      gcpBatchCallPaths.batchMonitoringLogPath, localMonitoringLogPath, workingDisk, optional = false, secondary = false,
      contentType = plainTextContentType)
  }


  // Primary entry point for cromwell to run GCP Batch job
  override def executeAsync(): Future[ExecutionHandle] = {

    def generateInputOutputParameters: Future[InputOutputParameters] = Future.fromTry(Try {
      val rcFileOutput = GcpBatchFileOutput(returnCodeFilename, returnCodeGcsPath, DefaultPathBuilder.get(returnCodeFilename), workingDisk, optional = false, secondary = false,
        contentType = plainTextContentType)

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
        GcpBatchFileOutput(s.name, returnCodeGcsPath.sibling(s.filename), DefaultPathBuilder.get(s.filename),
          workingDisk, optional = false, secondary = false, uploadPeriod = batchAttributes.logFlushPeriod, contentType = plainTextContentType)
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


    println(jobDescriptor.fullyQualifiedInputs)
    val file = jobDescriptor.localInputs
    println(file.get("test"))
    val gcpBatchParameters = CreateGcpBatchParameters(jobDescriptor = jobDescriptor, runtimeAttributes = runtimeAttributes, batchAttributes = batchAttributes, dockerImage = jobDockerImage, projectId = batchAttributes.project, region = batchAttributes.location)

    val runBatchResponse = for {
      //_ <- evaluateRuntimeAttributes
      _ <- uploadScriptFile()
      customLabels <- Future.fromTry(GoogleLabels.fromWorkflowOptions(workflowDescriptor.workflowOptions))
      jesParameters <- generateInputOutputParameters
      createParameters = createPipelineParameters(jesParameters)
      drsLocalizationManifestCloudPath = jobPaths.callExecutionRoot / GcpBatchJobPaths.DrsLocalizationManifestName
      // _ <- uploadDrsLocalizationManifest(createParameters, drsLocalizationManifestCloudPath)
      gcsTransferConfiguration = initializationData.gcpBatchConfiguration.batchAttributes.gcsTransferConfiguration
      gcsTransferLibraryCloudPath = jobPaths.callExecutionRoot / GcpBatchJobPaths.GcsTransferLibraryName
      transferLibraryContainerPath = createParameters.commandScriptContainerPath.sibling(GcsTransferLibraryName)
      _ <- uploadGcsTransferLibrary(createParameters, gcsTransferLibraryCloudPath, gcsTransferConfiguration)
      gcsLocalizationScriptCloudPath = jobPaths.callExecutionRoot / GcpBatchJobPaths.GcsLocalizationScriptName
      gcsDelocalizationScriptCloudPath = jobPaths.callExecutionRoot / GcpBatchJobPaths.GcsDelocalizationScriptName
      transferLibraryContainerPath = createParameters.commandScriptContainerPath.sibling(GcsTransferLibraryName)
      _ <- uploadGcsTransferLibrary(createParameters, gcsTransferLibraryCloudPath, gcsTransferConfiguration)
      _ = backendSingletonActor ! GcpBatchRequest(workflowId, jobName = jobTemp, gcpBatchCommand, gcpBatchParameters)
      runId = StandardAsyncJob(jobTemp)

    }
    yield runId

    runBatchResponse map { runId => PendingExecutionHandle(jobDescriptor, runId, Option(Run(runId)), previousState = None) }

  }

  override def reconnectAsync(jobId: StandardAsyncJob): Future[ExecutionHandle] = {
    log.info("reconnect async runs") // in for debugging remove later
    val handle = PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState](jobDescriptor, jobId, Option(Run(jobId)), previousState = None)
    Future.successful(handle)
  }

  override lazy val pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(5
    .second, 5
    .minutes, 1.1)

  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(
    initialInterval = 5
      .seconds, maxInterval = 20
      .seconds, multiplier = 1.1)

  protected def sendIncrementMetricsForReferenceFiles(referenceInputFilesOpt: Option[Set[GcpBatchInput]]): Unit = {
    referenceInputFilesOpt match {
      case Some(referenceInputFiles) =>
        referenceInputFiles.foreach { referenceInputFile =>
          increment(NonEmptyList.of("referencefiles", referenceInputFile.relativeHostPath.pathAsString))
        }
      case _ =>
      // do nothing - reference disks feature is either not configured in Cromwell or disabled in workflow options
    }
  }

  protected def sendIncrementMetricsForDockerImageCache(dockerImageCacheDiskOpt: Option[String],
                                                        dockerImageAsSpecifiedByUser: String,
                                                        isDockerImageCacheUsageRequested: Boolean): Unit = {
    (isDockerImageCacheUsageRequested, dockerImageCacheDiskOpt) match {
      case (true, None) => increment(NonEmptyList("docker", List("image", "cache", "image_not_in_cache", dockerImageAsSpecifiedByUser)))
      case (true, Some(_)) => increment(NonEmptyList("docker", List("image", "cache", "used_image_from_cache", dockerImageAsSpecifiedByUser)))
      case (false, Some(_)) => increment(NonEmptyList("docker", List("image", "cache", "cached_image_not_used", dockerImageAsSpecifiedByUser)))
      case _ => // docker image cache not requested and image is not in cache anyway - do nothing
    }
  }

  override def pollStatusAsync(handle: GcpBatchPendingExecutionHandle): Future[RunStatus] = {

    val jobId = handle.pendingJob.jobId

    val job = handle.runInfo match {
      case Some(actualJob) => actualJob
      case None =>
        throw new RuntimeException(
          s"pollStatusAsync called but job not available. This should not happen. Job Id $jobId"
        )
    }

    log.info(s"started polling for $job with jobId $jobId")

    //super[GcpBatchStatusRequestClient].pollStatus(workflowId, handle.pendingJob, jobTemp)

    //temporary. Added to resolve issue with poll async starter before job submitted.
    implicit val timeout: Timeout = Timeout(60.seconds) //had to set to high amount for some reason.  Otherwise would not finish with low value
    val futureResult = backendSingletonActor ? BatchJobAsk(jobId)
    val result = Await.result(futureResult, timeout.duration).asInstanceOf[String]
    log.info(result)

    try {
      val gcpBatchPoll = new GcpBatchJobGetRequest
      val result = gcpBatchPoll.GetJob(jobId, batchAttributes.project, batchAttributes.location)
      RunStatus.fromJobStatus(result)

    }
    catch {
      case nfe: NotFoundException => //added to account for job not found errors because polling async happens before job is submitted
        nfe.printStackTrace()
        Future.successful(Running)
      case ee: ExecutionException => //added to account for job not found errors because polling async happens before job is submitted
        ee.printStackTrace()
        Future.successful(Running)
    }

  }

  override def isTerminal(runStatus: RunStatus): Boolean = {
    //runStatus.isTerminal

    runStatus match {
      case jobSucceeded: RunStatus.Succeeded =>
        log.info("isTerminal match Succeeded running with status {}", jobSucceeded)
        true
      case jobFailed: RunStatus.Failed =>
        log.info("isTerminal match Failed with status {}", jobFailed)
        true
      case _: TerminalRunStatus =>
        val tempTermStatus = runStatus.toString
        log.info(f"isTerminal match TerminalRunStatus running with status $tempTermStatus")
        true
      case other =>
        log.info(f"isTerminal match _ running with status $other")
        false
    }
  }

  override def isDone(runStatus: RunStatus): Boolean = {
    runStatus match {
      case _: RunStatus.Succeeded =>
        log.info("GCP batch job succeeded matched isDone")
        true
      case _: RunStatus.Failed =>
        log.info("GCP Job failed and matched isDone")
        true
      case _ =>
        log.info("did not match isDone")
        false //throw new RuntimeException(s"Cromwell programmer blunder: isSuccess was called on an incomplete RunStatus ($runStatus).")
    }
  }

  override def getTerminalEvents(runStatus: RunStatus): Seq[ExecutionEvent] = {
    runStatus match {
      case successStatus: Succeeded => successStatus
        .eventList
      case unknown =>
        throw new RuntimeException(s"handleExecutionSuccess not called with RunStatus.Success. Instead got $unknown")
    }
  }

  override def getTerminalMetadata(runStatus: RunStatus): Map[String, Any] = {
    runStatus match {
      case _: TerminalRunStatus => Map()
      case unknown => throw new RuntimeException(s"Attempt to get terminal metadata from non terminal status: $unknown")
    }
  }

  override val gcpBatchActor: ActorRef = backendSingletonActor

  /*
  override def globParentDirectory(womGlobFile: WomGlobFile): Path = {
    val (_, disk) = relativePathAndAttachedDisk(womGlobFile.value, runtimeAttributes.disks)
    disk.mountPoint
  }
   */

  protected def googleProject(descriptor: BackendWorkflowDescriptor): String = {
    descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleProject, batchAttributes.project)
  }

  protected def computeServiceAccount(descriptor: BackendWorkflowDescriptor): String = {
    descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleComputeServiceAccount, batchAttributes.computeServiceAccount)
  }

  protected def fuseEnabled(descriptor: BackendWorkflowDescriptor): Boolean = {
    descriptor.workflowOptions.getBoolean(WorkflowOptionKeys.EnableFuse).toOption.getOrElse(batchAttributes.enableFuse)
  }

  protected def useDockerImageCache(descriptor: BackendWorkflowDescriptor): Boolean = {
    descriptor.workflowOptions.getBoolean(WorkflowOptionKeys.UseDockerImageCache).getOrElse(false)
  }

  override def cloudResolveWomFile(womFile: WomFile): WomFile = {
    womFile.mapFile { value =>
      getPath(value) match {
        case Success(drsPath: DrsPath) => DrsResolver.getSimpleGsUri(drsPath).unsafeRunSync().getOrElse(value)
        case Success(path) => path.pathAsString
        case _ => value
      }
    }
  }

  override def mapCommandLineWomFile(womFile: WomFile): WomFile = {
    womFile.mapFile { value =>
      (getPath(value), asAdHocFile(womFile)) match {
        case (Success(gcsPath: GcsPath), Some(adHocFile)) =>
          // Ad hoc files will be placed directly at the root ("/cromwell_root/ad_hoc_file.txt") unlike other input files
          // for which the full path is being propagated ("/cromwell_root/path/to/input_file.txt")
          workingDisk.mountPoint.resolve(adHocFile.alternativeName.getOrElse(gcsPath.name)).pathAsString
        case (Success(path@(_: GcsPath | _: HttpPath)), _) =>
          workingDisk.mountPoint.resolve(path.pathWithoutScheme).pathAsString
        case (Success(drsPath: DrsPath), _) =>
          val filePath = DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
          workingDisk.mountPoint.resolve(filePath).pathAsString
        case (Success(sraPath: SraPath), _) =>
          workingDisk.mountPoint.resolve(s"sra-${sraPath.accession}/${sraPath.pathWithoutScheme}").pathAsString
        case _ => value
      }
    }
  }

  override def mapCommandLineJobInputWomFile(womFile: WomFile): WomFile = {
    womFile.mapFile(value =>
      getPath(value) match {
        case Success(gcsPath: GcsPath) => workingDisk.mountPoint.resolve(gcsPath.pathWithoutScheme).pathAsString
        case Success(drsPath: DrsPath) =>
          val filePath = DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
          workingDisk.mountPoint.resolve(filePath).pathAsString
        case _ => value
      }
    )
  }






}




