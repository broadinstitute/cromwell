package cromwell.backend.impl.jes

import java.net.SocketTimeoutException
import java.nio.file.{Path, Paths}

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.event.LoggingReceive
import better.files._
import cats.instances.future._
import cats.syntax.functor._
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import cromwell.backend.BackendJobExecutionActor.{AbortedResponse, BackendJobExecutionResponse}
import cromwell.backend.BackendLifecycleActor.AbortJobCommand
import cromwell.backend.async.AsyncBackendJobExecutionActor.{ExecutionMode, JobId}
import cromwell.backend.async.{AbortedExecutionHandle, AsyncBackendJobExecutionActor, ExecutionHandle, FailedNonRetryableExecutionHandle, FailedRetryableExecutionHandle, NonRetryableExecution, SuccessfulExecutionHandle}
import cromwell.backend.impl.jes.JesJobExecutionActor.JesOperationIdKey
import cromwell.backend.impl.jes.RunStatus.TerminalRunStatus
import cromwell.backend.impl.jes.io._
import cromwell.backend.impl.jes.statuspolling.JesPollingActorClient
import cromwell.backend.wdl.OutputEvaluator
import cromwell.backend.{BackendJobDescriptor, BackendWorkflowDescriptor, PreemptedException}
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core._
import cromwell.core.logging.JobLogging
import cromwell.core.path.proxy.PathProxy
import cromwell.core.retry.{Retry, SimpleExponentialBackoff}
import cromwell.core.path.PathImplicits._
import cromwell.services.keyvalue.KeyValueServiceActor._
import cromwell.services.metadata._
import wdl4s._
import wdl4s.expression.NoFunctions
import wdl4s.values._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object JesAsyncBackendJobExecutionActor {

  def props(jobDescriptor: BackendJobDescriptor,
            completionPromise: Promise[BackendJobExecutionResponse],
            jesWorkflowInfo: JesConfiguration,
            initializationData: JesBackendInitializationData,
            serviceRegistryActor: ActorRef,
            jesBackendSingletonActor: ActorRef): Props = {
    Props(new JesAsyncBackendJobExecutionActor(jobDescriptor,
      completionPromise,
      jesWorkflowInfo,
      initializationData,
      serviceRegistryActor,
      jesBackendSingletonActor)).withDispatcher(BackendDispatcher)
  }

  object WorkflowOptionKeys {
    val MonitoringScript = "monitoring_script"
    val GoogleProject = "google_project"
    val GoogleComputeServiceAccount = "google_compute_service_account"
  }


  private val ExtraConfigParamName = "__extra_config_gcs_path"

  /**
    * Representing a running JES execution, instances of this class are never Done and it is never okay to
    * ask them for results.
    */
  case class JesPendingExecutionHandle(jobDescriptor: BackendJobDescriptor,
                                       jesOutputs: Set[JesFileOutput],
                                       run: Run,
                                       previousStatus: Option[RunStatus]) extends ExecutionHandle {
    override val isDone = false
    override val result = NonRetryableExecution(new IllegalStateException("JesPendingExecutionHandle cannot yield a result"))
  }

  case class JesJobId(operationId: String) extends JobId
}


class JesAsyncBackendJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                                       override val completionPromise: Promise[BackendJobExecutionResponse],
                                       override val jesConfiguration: JesConfiguration,
                                       override val initializationData: JesBackendInitializationData,
                                       override val serviceRegistryActor: ActorRef,
                                       val jesBackendSingletonActor: ActorRef)
  extends Actor with ActorLogging with AsyncBackendJobExecutionActor with JesJobCachingActorHelper with JobLogging with JesPollingActorClient {

  import JesAsyncBackendJobExecutionActor._

  override val pollingActor = jesBackendSingletonActor

  override lazy val pollBackOff = SimpleExponentialBackoff(
    initialInterval = 30 seconds, maxInterval = jesAttributes.maxPollingInterval seconds, multiplier = 1.1)

  override lazy val executeOrRecoverBackOff = SimpleExponentialBackoff(
    initialInterval = 3 seconds, maxInterval = 20 seconds, multiplier = 1.1)

  override lazy val workflowDescriptor = jobDescriptor.workflowDescriptor

  private lazy val call = jobDescriptor.key.call

  override lazy val retryable = jobDescriptor.key.attempt <= runtimeAttributes.preemptible
  private lazy val cmdInput =
    JesFileInput(ExecParamName, jesCallPaths.script.toRealString, Paths.get(jesCallPaths.scriptFilename), workingDisk)
  private lazy val jesCommandLine = s"/bin/bash ${cmdInput.containerPath}"
  private lazy val rcJesOutput = JesFileOutput(returnCodeFilename, returnCodeGcsPath.toRealString, Paths.get(returnCodeFilename), workingDisk)

  private lazy val standardParameters = Seq(rcJesOutput)
  private lazy val returnCodeContents = Try(File(returnCodeGcsPath).contentAsString)
  private lazy val dockerConfiguration = jesConfiguration.dockerCredentials
  private lazy val tag = s"${this.getClass.getSimpleName} [UUID(${workflowId.shortString}):${jobDescriptor.key.tag}]"

  private var runId: Option[String] = None

  def jesReceiveBehavior: Receive = LoggingReceive {
    case AbortJobCommand =>
      runId foreach { id =>
        Try(Run(id, initializationData.genomics).abort()) match {
          case Success(_) => jobLogger.info("{} Aborted {}", tag: Any, id)
          case Failure(ex) => jobLogger.warn("{} Failed to abort {}: {}", tag, id, ex.getMessage)
        }
      }
      context.parent ! AbortedResponse(jobDescriptor.key)
      context.stop(self)
    case KvPutSuccess(_) => // expected after the KvPut for the operation ID
  }

  override def receive: Receive = pollingActorClientReceive orElse jesReceiveBehavior orElse super.receive

  private def gcsAuthParameter: Option[JesInput] = {
    if (jesAttributes.auths.gcs.requiresAuthFile || dockerConfiguration.isDefined)
      Option(JesLiteralInput(ExtraConfigParamName, jesCallPaths.gcsAuthFilePath.toRealString))
    else None
  }

  private lazy val callContext = CallContext(
    callRootPath,
    jesStdoutFile.toRealString,
    jesStderrFile.toRealString
  )

  private[jes] lazy val callEngineFunctions = new JesExpressionFunctions(List(jesCallPaths.gcsPathBuilder), callContext)

  /**
    * Takes two arrays of remote and local WDL File paths and generates the necessary JesInputs.
    */
  private def jesInputsFromWdlFiles(jesNamePrefix: String,
                                    remotePathArray: Seq[WdlFile],
                                    localPathArray: Seq[WdlFile],
                                    jobDescriptor: BackendJobDescriptor): Iterable[JesInput] = {
    (remotePathArray zip localPathArray zipWithIndex) flatMap {
      case ((remotePath, localPath), index) =>
        Seq(JesFileInput(s"$jesNamePrefix-$index", remotePath.valueString, Paths.get(localPath.valueString), workingDisk))
    }
  }

  /**
    * Turns WdlFiles into relative paths.  These paths are relative to the working disk
    *
    * relativeLocalizationPath("foo/bar.txt") -> "foo/bar.txt"
    * relativeLocalizationPath("gs://some/bucket/foo.txt") -> "some/bucket/foo.txt"
    */
  private def relativeLocalizationPath(file: WdlFile): WdlFile = {
    getPath(file.value) match {
      case Success(path) =>
        val value: WdlSource = path.toUri.getHost + path.toUri.getPath
        WdlFile(value, file.isGlob)
      case _ => file
    }
  }

  private[jes] def generateJesInputs(jobDescriptor: BackendJobDescriptor): Set[JesInput] = {

    val writeFunctionFiles = call.task.evaluateFilesFromCommand(jobDescriptor.fullyQualifiedInputs, callEngineFunctions) map {
      case (expression, file) =>  expression.toWdlString.md5SumShort -> Seq(file)
    }

    /* Collect all WdlFiles from inputs to the call */
    val callInputFiles: Map[FullyQualifiedName, Seq[WdlFile]] = jobDescriptor.fullyQualifiedInputs mapValues { _.collectAsSeq { case w: WdlFile => w } }

    val inputs = (callInputFiles ++ writeFunctionFiles) flatMap {
      case (name, files) => jesInputsFromWdlFiles(name, files, files.map(relativeLocalizationPath), jobDescriptor)
    }

    inputs.toSet
  }

  /**
    * Given a path (relative or absolute), returns a (Path, JesAttachedDisk) tuple where the Path is
    * relative to the AttachedDisk's mount point
    *
    * @throws Exception if the `path` does not live in one of the supplied `disks`
    */
  private def relativePathAndAttachedDisk(path: String, disks: Seq[JesAttachedDisk]): (Path, JesAttachedDisk) = {
    val absolutePath = Paths.get(path) match {
      case p if !p.isAbsolute => JesWorkingDisk.MountPoint.resolve(p)
      case p => p
    }

    disks.find(d => absolutePath.startsWith(d.mountPoint)) match {
      case Some(disk) => (disk.mountPoint.relativize(absolutePath), disk)
      case None =>
        throw new Exception(s"Absolute path $path doesn't appear to be under any mount points: ${disks.map(_.toString).mkString(", ")}")
    }
  }

  /**
    * If the desired reference name is too long, we don't want to break JES or risk collisions by arbitrary truncation. So,
    * just use a hash. We only do this when needed to give better traceability in the normal case.
    */
  private def makeSafeJesReferenceName(referenceName: String) = {
    if (referenceName.length <= 127) referenceName else referenceName.md5Sum
  }

  private[jes] def findGlobOutputs(jobDescriptor: BackendJobDescriptor): Set[WdlGlobFile] = {
    val globOutputs = (call.task.findOutputFiles(jobDescriptor.fullyQualifiedInputs, NoFunctions) map relativeLocalizationPath) collect {
      case glob: WdlGlobFile => glob
    }
    globOutputs.distinct.toSet
  }

  private[jes] def generateJesOutputs(jobDescriptor: BackendJobDescriptor): Set[JesFileOutput] = {
    val wdlFileOutputs = call.task.findOutputFiles(jobDescriptor.fullyQualifiedInputs, NoFunctions) map relativeLocalizationPath

    val outputs = wdlFileOutputs.distinct flatMap { wdlFile =>
      wdlFile match {
        case singleFile: WdlSingleFile => List(generateJesSingleFileOutputs(singleFile))
        case globFile: WdlGlobFile => generateJesGlobFileOutputs(globFile)
      }
    }

    outputs.toSet
  }

  private def generateJesSingleFileOutputs(wdlFile: WdlSingleFile): JesFileOutput = {
    val destination = callRootPath.resolve(wdlFile.value.stripPrefix("/")).toRealString
    val (relpath, disk) = relativePathAndAttachedDisk(wdlFile.value, runtimeAttributes.disks)
    JesFileOutput(makeSafeJesReferenceName(wdlFile.value), destination, relpath, disk)
  }

  private def generateJesGlobFileOutputs(wdlFile: WdlGlobFile): List[JesFileOutput] = {
    val globName = callEngineFunctions.globName(wdlFile.value)
    val globDirectory = globName + "/"
    val globListFile = globName + ".list"
    val gcsGlobDirectoryDestinationPath = callRootPath.resolve(globDirectory).toRealString
    val gcsGlobListFileDestinationPath = callRootPath.resolve(globListFile).toRealString

    val (_, globDirectoryDisk) = relativePathAndAttachedDisk(wdlFile.value, runtimeAttributes.disks)

    // We need both the glob directory and the glob list:
    List(
      // The glob directory:
      JesFileOutput(makeSafeJesReferenceName(globDirectory), gcsGlobDirectoryDestinationPath, Paths.get(globDirectory + "*"), globDirectoryDisk),
      // The glob list file:
      JesFileOutput(makeSafeJesReferenceName(globListFile), gcsGlobListFileDestinationPath, Paths.get(globListFile), globDirectoryDisk)
    )
  }

  private def instantiateCommand: Try[String] = {
    val backendInputs = jobDescriptor.inputDeclarations mapValues gcsPathToLocal
    jobDescriptor.call.task.instantiateCommand(backendInputs, callEngineFunctions, valueMapper = gcsPathToLocal)
  }

  private def uploadCommandScript(command: String, withMonitoring: Boolean, globFiles: Set[WdlGlobFile]): Future[Unit] = {
    val monitoring = if (withMonitoring) {
      s"""|touch $JesMonitoringLogFile
          |chmod u+x $JesMonitoringScript
          |$JesMonitoringScript > $JesMonitoringLogFile &""".stripMargin
    } else ""

    val tmpDir = File(JesWorkingDisk.MountPoint)./("tmp").path
    val rcPath = File(JesWorkingDisk.MountPoint)./(returnCodeFilename).path

    def globManipulation(globFile: WdlGlobFile) = {

      val globDir = callEngineFunctions.globName(globFile.value)
      val (_, disk) = relativePathAndAttachedDisk(globFile.value, runtimeAttributes.disks)
      val globDirectory = Paths.get(s"${disk.mountPoint.toAbsolutePath}/$globDir/")
      val globList = Paths.get(s"${disk.mountPoint.toAbsolutePath}/$globDir.list")

      s"""
        |mkdir $globDirectory
        |ln ${globFile.value} $globDirectory
        |ls -1 $globDirectory > $globList
      """.stripMargin
    }

    val globManipulations = globFiles.map(globManipulation).mkString("\n")

    val fileContent =
      s"""
         |#!/bin/bash
         |export _JAVA_OPTIONS=-Djava.io.tmpdir=$tmpDir
         |export TMPDIR=$tmpDir
         |$monitoring
         |(
         |cd ${JesWorkingDisk.MountPoint}
         |$command
         |$globManipulations
         |)
         |echo $$? > $rcPath
       """.stripMargin.trim

    Future(File(jesCallPaths.script).write(fileContent)) void
  }

  private def googleProject(descriptor: BackendWorkflowDescriptor): String = {
    descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleProject, jesAttributes.project)
  }

  private def computeServiceAccount(descriptor: BackendWorkflowDescriptor): String = {
    descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleComputeServiceAccount, jesAttributes.computeServiceAccount)
  }

  private def createJesRun(jesParameters: Seq[JesParameter], runIdForResumption: Option[String]): Future[Run] = {

    def createRun() = Future(Run(
      runIdForResumption,
      jobDescriptor = jobDescriptor,
      runtimeAttributes = runtimeAttributes,
      callRootPath = callRootPath.toRealString,
      commandLine = jesCommandLine,
      logFileName = jesLogFilename,
      jesParameters,
      googleProject(jobDescriptor.workflowDescriptor),
      computeServiceAccount(jobDescriptor.workflowDescriptor),
      retryable,
      initializationData.genomics
    ))

    implicit val system = context.system
    Retry.withRetry(
      createRun,
      isTransient = isTransientJesException,
      isFatal = isFatalJesException
    ) andThen {
      case Success(run) =>
        // If this execution represents a resumption don't publish the operation ID since clearly it is already persisted.
        runId = Option(run.runId)
        if (runIdForResumption.isEmpty) {
          serviceRegistryActor ! KvPut(KvPair(ScopedKey(jobDescriptor.workflowDescriptor.id,
            KvJobKey(jobDescriptor.key.call.fullyQualifiedName, jobDescriptor.key.index, jobDescriptor.key.attempt),
            JesOperationIdKey), runId))
        }
    }
  }

  protected def runWithJes(command: String,
                           jesInputs: Set[JesInput],
                           jesOutputs: Set[JesFileOutput],
                           runIdForResumption: Option[String],
                           withMonitoring: Boolean): Future[ExecutionHandle] = {

    tellStartMetadata()

    val jesParameters = standardParameters ++ gcsAuthParameter ++ jesInputs ++ jesOutputs

    val jesJobSetup = for {
      _ <- uploadCommandScript(command, withMonitoring, findGlobOutputs(jobDescriptor))
      run <- createJesRun(jesParameters, runIdForResumption)
      _ = tellMetadata(Map(CallMetadataKeys.JobId -> run.runId))
    } yield run

    jesJobSetup map { run => JesPendingExecutionHandle(jobDescriptor, jesOutputs, run, previousStatus = None) }
  }

  override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    // Force runtimeAttributes to evaluate so we can fail quickly now if we need to:
    Try(runtimeAttributes) match {
      case Success(_) => startExecuting(monitoringOutput, mode)
      case Failure(e) => Future.successful(FailedNonRetryableExecutionHandle(e, None))
    }
  }

  private def startExecuting(monitoringOutput: Option[JesFileOutput], mode: ExecutionMode): Future[ExecutionHandle] = {
    val jesInputs: Set[JesInput] = generateJesInputs(jobDescriptor) ++ monitoringScript + cmdInput
    val jesOutputs: Set[JesFileOutput] = generateJesOutputs(jobDescriptor) ++ monitoringOutput

    instantiateCommand match {
      case Success(command) => runWithJes(command, jesInputs, jesOutputs, mode.jobId.collectFirst { case j: JesJobId => j.operationId }, monitoringScript.isDefined)
      case Failure(ex: SocketTimeoutException) => Future.successful(FailedNonRetryableExecutionHandle(ex))
      case Failure(ex) => Future.successful(FailedNonRetryableExecutionHandle(ex))
    }
  }

  /**
    * Update the ExecutionHandle
    */
  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    previous match {
      case handle: JesPendingExecutionHandle =>
        jobLogger.debug(s"$tag Polling JES Job ${handle.run.runId}")
        pollStatus(handle.run) map updateExecutionHandleSuccess(handle) recover updateExecutionHandleFailure(handle) flatten
      case f: FailedNonRetryableExecutionHandle => f.future
      case s: SuccessfulExecutionHandle => s.future
      case badHandle => Future.failed(new IllegalArgumentException(s"Unexpected execution handle: $badHandle"))
    }
  }

  private def updateExecutionHandleFailure(oldHandle: JesPendingExecutionHandle): PartialFunction[Throwable, Future[ExecutionHandle]] = {
    case e: GoogleJsonResponseException if e.getStatusCode == 404 =>
      jobLogger.error(s"$tag JES Job ID ${oldHandle.run.runId} has not been found, failing call")
      FailedNonRetryableExecutionHandle(e).future
    case e: Exception =>
      // Log exceptions and return the original handle to try again.
      jobLogger.warn(s"Caught exception, retrying", e)
      oldHandle.future
    case e: Error => Future.failed(e) // JVM-ending calamity.
    case throwable =>
      // Someone has subclassed Throwable directly?
      FailedNonRetryableExecutionHandle(throwable).future
  }

  private def updateExecutionHandleSuccess(oldHandle: JesPendingExecutionHandle)(status: RunStatus): Future[ExecutionHandle] = {
    val previousStatus = oldHandle.previousStatus
    if (!(previousStatus contains status)) {
      // If this is the first time checking the status, we log the transition as '-' to 'currentStatus'. Otherwise
      // just use the state names.
      val prevStateName = previousStatus map { _.toString } getOrElse "-"
      jobLogger.info(s"$tag Status change from $prevStateName to $status")
      tellMetadata(Map("backendStatus" -> status))
    }

    status match {
      case s: TerminalRunStatus =>
        val metadata = Map(
          JesMetadataKeys.MachineType -> s.machineType.getOrElse("unknown"),
          JesMetadataKeys.InstanceName -> s.instanceName.getOrElse("unknown"),
          JesMetadataKeys.Zone -> s.zone.getOrElse("unknown")
        )

        tellMetadata(metadata)
        executionResult(s, oldHandle)
      case s => oldHandle.copy(previousStatus = Option(s)).future // Copy the current handle with updated previous status.

    }
  }

  /**
    * Fire and forget start info to the metadata service
    */
  private def tellStartMetadata() = tellMetadata(metadataKeyValues)

  /**
    * Fire and forget info to the metadata service
    */
  def tellMetadata(metadataKeyValues: Map[String, Any]): Unit = {
    import cromwell.services.metadata.MetadataService.implicits.MetadataAutoPutter
    serviceRegistryActor.putMetadata(jobDescriptor.workflowDescriptor.id, Option(jobDescriptor.key), metadataKeyValues)
  }

  private[jes] def wdlValueToGcsPath(jesOutputs: Set[JesFileOutput])(value: WdlValue): WdlValue = {
    def toGcsPath(wdlFile: WdlFile) = jesOutputs collectFirst {
      case o if o.name == makeSafeJesReferenceName(wdlFile.valueString) => WdlFile(o.gcs)
    } getOrElse value

    value match {
      case wdlArray: WdlArray => wdlArray map wdlValueToGcsPath(jesOutputs)
      case wdlMap: WdlMap => wdlMap map {
        case (k, v) => wdlValueToGcsPath(jesOutputs)(k) -> wdlValueToGcsPath(jesOutputs)(v)
      }
      case file: WdlFile => toGcsPath(file)
      case other => other
    }
  }

  private def postProcess: Try[CallOutputs] = {
    def wdlValueToSuccess(value: WdlValue): Try[WdlValue] = Success(value)

    OutputEvaluator.evaluateOutputs(
      jobDescriptor,
      callEngineFunctions,
      (wdlValueToSuccess _).compose(wdlValueToGcsPath(generateJesOutputs(jobDescriptor)))
    )
  }

  private def handleSuccess(outputMappings: Try[CallOutputs], returnCode: Int, jobDetritusFiles: Map[String, Path], executionHandle: ExecutionHandle, events: Seq[ExecutionEvent]): ExecutionHandle = {
    outputMappings match {
      case Success(outputs) => SuccessfulExecutionHandle(outputs, returnCode, jobDetritusFiles, events)
      case Failure(ex: CromwellAggregatedException) if ex.throwables collectFirst { case s: SocketTimeoutException => s } isDefined =>
        // Return the execution handle in this case to retry the operation
        executionHandle
      case Failure(ex) => FailedNonRetryableExecutionHandle(ex)
    }
  }

  private def extractErrorCodeFromErrorMessage(errorMessage: String): Int = {
    errorMessage.substring(0, errorMessage.indexOf(':')).toInt
  }

  private def preempted(errorCode: Int, errorMessage: List[String]): Boolean = {
    def isPreemptionCode(code: Int) = code == 13 || code == 14

    try {
      errorCode == 10 && errorMessage.exists(e => isPreemptionCode(extractErrorCodeFromErrorMessage(e))) && preemptible
    } catch {
      case _: NumberFormatException | _: StringIndexOutOfBoundsException =>
        jobLogger.warn(s"Unable to parse JES error code from error messages: [{}], assuming this was not a preempted VM.", errorMessage.mkString(", "))
        false
    }
  }

  private def handleFailure(errorCode: Int, errorMessage: List[String]) = {
    import lenthall.numeric.IntegerUtil._

    val taskName = s"${workflowDescriptor.id}:${call.unqualifiedName}"
    val attempt = jobDescriptor.key.attempt

    if (errorMessage.exists(_.contains("Operation canceled at")))  {
      AbortedExecutionHandle.future
    } else if (preempted(errorCode, errorMessage)) {
      val preemptedMsg = s"Task $taskName was preempted for the ${attempt.toOrdinal} time."

      if (attempt < maxPreemption) {
        val e = PreemptedException(
          s"""$preemptedMsg The call will be restarted with another preemptible VM (max preemptible attempts number is $maxPreemption).
             |Error code $errorCode. Message: $errorMessage""".stripMargin
        )
        FailedRetryableExecutionHandle(e, None).future
      } else {
        val e = PreemptedException(
          s"""$preemptedMsg The maximum number of preemptible attempts ($maxPreemption) has been reached. The call will be restarted with a non-preemptible VM.
             |Error code $errorCode. Message: $errorMessage)""".stripMargin)
        FailedRetryableExecutionHandle(e, None).future
      }
    } else {
      val id = workflowDescriptor.id
      val name = jobDescriptor.call.unqualifiedName
      val message = if (errorMessage.isEmpty) "null" else errorMessage.mkString(", ")
      val exception = new RuntimeException(s"Task $id:$name failed: error code $errorCode. Message: $message")
      FailedNonRetryableExecutionHandle(exception, None).future
    }
  }

  private[jes] def executionResult(status: TerminalRunStatus, handle: JesPendingExecutionHandle)
                                  (implicit ec: ExecutionContext): Future[ExecutionHandle] = Future {
    try {
      lazy val stderrLength: Long = File(jesStderrFile).size
      lazy val returnCode = returnCodeContents map { _.trim.toInt }
      lazy val continueOnReturnCode = runtimeAttributes.continueOnReturnCode

      status match {
        case _: RunStatus.Success if runtimeAttributes.failOnStderr && stderrLength.intValue > 0 =>
          // returnCode will be None if it couldn't be downloaded/parsed, which will yield a null in the DB
          FailedNonRetryableExecutionHandle(new RuntimeException(
            s"execution failed: stderr has length $stderrLength"), returnCode.toOption).future
        case _: RunStatus.Success if returnCodeContents.isFailure =>
          val exception = returnCode.failed.get
          jobLogger.warn(s"could not download return code file, retrying", exception)
          // Return handle to try again.
          handle.future
        case _: RunStatus.Success if returnCode.isFailure =>
          FailedNonRetryableExecutionHandle(new RuntimeException(
            s"execution failed: could not parse return code as integer: ${returnCodeContents.get}")).future
        case _: RunStatus.Success if !continueOnReturnCode.continueFor(returnCode.get) =>
          val badReturnCodeMessage = s"Call ${jobDescriptor.key.tag}: return code was ${returnCode.getOrElse("(none)")}"
          FailedNonRetryableExecutionHandle(new RuntimeException(badReturnCodeMessage), returnCode.toOption).future
        case success: RunStatus.Success =>
          handleSuccess(postProcess, returnCode.get, jesCallPaths.detritusPaths, handle, success.eventList).future
        case RunStatus.Failed(errorCode, errorMessage, _, _, _, _) => handleFailure(errorCode, errorMessage)
      }
    } catch {
      case e: Exception =>
        jobLogger.warn("Caught exception trying to download result, retrying", e)
        // Return the original handle to try again.
        handle.future
    }
  } flatten

  /**
    * Takes a path in GCS and comes up with a local path which is unique for the given GCS path.
    *
    * Matches the path generated via relativeLocalizationPath and passed in as JesFileInput.local.
    *
    * @param mountPoint The mount point for inputs
    * @param gcsPath The input path
    * @return A path which is unique per input path
    */
  private def localFilePathFromCloudStoragePath(mountPoint: Path, gcsPath: CloudStoragePath): Path = {
    mountPoint.resolve(gcsPath.bucket()).resolve(gcsPath.toUri.getPath.stripPrefix("/"))
  }

  /**
    * Takes a single WdlValue and maps google cloud storage (GCS) paths into an appropriate local file path.
    * If the input is not a WdlFile, or the WdlFile is not a GCS path, the mapping is a noop.
    *
    * @param wdlValue the value of the input
    * @return a new FQN to WdlValue pair, with WdlFile paths modified if appropriate.
    */
  private[jes] def gcsPathToLocal(wdlValue: WdlValue): WdlValue = {
    wdlValue match {
      case wdlFile: WdlFile =>
        getPath(wdlFile.valueString) match {
          case Success(gcsPath: CloudStoragePath) =>
            WdlFile(localFilePathFromCloudStoragePath(workingDisk.mountPoint, gcsPath).toString, wdlFile.isGlob)
          case Success(proxy: PathProxy) =>
            proxy.unbox(classOf[CloudStoragePath]) map { gcsPath =>
              WdlFile(localFilePathFromCloudStoragePath(workingDisk.mountPoint, gcsPath).toString, wdlFile.isGlob)
            } getOrElse wdlValue
          case _ => wdlValue
        }
      case wdlArray: WdlArray => wdlArray map gcsPathToLocal
      case wdlMap: WdlMap => wdlMap map { case (k, v) => gcsPathToLocal(k) -> gcsPathToLocal(v) }
      case _ => wdlValue
    }
  }

  protected implicit def ec: ExecutionContext = context.dispatcher
}
