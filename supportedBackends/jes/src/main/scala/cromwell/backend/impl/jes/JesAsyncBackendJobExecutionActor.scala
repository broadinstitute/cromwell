package cromwell.backend.impl.jes

import java.net.SocketTimeoutException
import java.nio.file.{Path, Paths}

import akka.actor.{Actor, ActorLogging, ActorRef}
import better.files._
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend._
import cromwell.backend.async.{AbortedExecutionHandle, ExecutionHandle, FailedNonRetryableExecutionHandle, FailedRetryableExecutionHandle, PendingExecutionHandle, SuccessfulExecutionHandle}
import cromwell.backend.impl.jes.RunStatus.TerminalRunStatus
import cromwell.backend.impl.jes.io._
import cromwell.backend.impl.jes.statuspolling.JesPollingActorClient
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.backend.validation.ContinueOnReturnCode
import cromwell.backend.wdl.OutputEvaluator
import cromwell.core._
import cromwell.core.logging.JobLogging
import cromwell.core.path.PathFactory._
import cromwell.core.path.PathImplicits._
import cromwell.core.path.proxy.PathProxy
import cromwell.core.retry.SimpleExponentialBackoff
import wdl4s._
import wdl4s.expression.{NoFunctions, WdlFunctions}
import wdl4s.values._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

case class JesAsyncExecutionActorParams
(
  override val jobDescriptor: BackendJobDescriptor,
  jesConfiguration: JesConfiguration,
  jesBackendInitializationData: JesBackendInitializationData,
  override val serviceRegistryActor: ActorRef,
  jesBackendSingletonActorOption: Option[ActorRef],
  override val completionPromise: Promise[BackendJobExecutionResponse]
) extends StandardAsyncExecutionActorParams {
  override val jobIdKey: String = JesJobExecutionActor.JesOperationIdKey
  override val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
  override val backendInitializationDataOption: Option[BackendInitializationData] = Option(jesBackendInitializationData)
}

object JesAsyncBackendJobExecutionActor {
  object WorkflowOptionKeys {
    val MonitoringScript = "monitoring_script"
    val GoogleProject = "google_project"
    val GoogleComputeServiceAccount = "google_compute_service_account"
  }

  type JesPendingExecutionHandle =
    PendingExecutionHandle[StandardAsyncJob, Run, RunStatus]

  private val ExtraConfigParamName = "__extra_config_gcs_path"
}

class JesAsyncBackendJobExecutionActor(val jesParams: JesAsyncExecutionActorParams)
  extends Actor with ActorLogging with BackendJobLifecycleActor with StandardAsyncExecutionActor
    with JesJobCachingActorHelper with JobLogging with JesPollingActorClient {

  import JesAsyncBackendJobExecutionActor._

  override val standardParams: StandardAsyncExecutionActorParams = jesParams
  override val jesConfiguration: JesConfiguration = jesParams.jesConfiguration
  override val initializationData: JesBackendInitializationData = {
    backendInitializationDataAs[JesBackendInitializationData]
  }

  val jesBackendSingletonActor: ActorRef =
    jesParams.jesBackendSingletonActorOption.getOrElse(
      throw new RuntimeException("JES Backend actor cannot exist without the JES backend singleton actor"))

  override type StandardAsyncRunInfo = Run

  override type StandardAsyncRunStatus = RunStatus

  override val pollingActor: ActorRef = jesBackendSingletonActor

  override lazy val pollBackOff = SimpleExponentialBackoff(
    initialInterval = 30 seconds, maxInterval = jesAttributes.maxPollingInterval seconds, multiplier = 1.1)

  override lazy val executeOrRecoverBackOff = SimpleExponentialBackoff(
    initialInterval = 3 seconds, maxInterval = 20 seconds, multiplier = 1.1)

  override lazy val workflowDescriptor: BackendWorkflowDescriptor = jobDescriptor.workflowDescriptor

  lazy val call: TaskCall = jobDescriptor.key.call

  override lazy val retryable: Boolean = jobDescriptor.key.attempt <= runtimeAttributes.preemptible
  private lazy val cmdInput =
    JesFileInput(ExecParamName, jesCallPaths.script.toRealString, Paths.get(jesCallPaths.scriptFilename), workingDisk)
  private lazy val jesCommandLine = s"/bin/bash ${cmdInput.containerPath}"
  private lazy val rcJesOutput = JesFileOutput(returnCodeFilename, returnCodeGcsPath.toRealString, Paths.get(returnCodeFilename), workingDisk)

  private lazy val standardParameters = Seq(rcJesOutput)

  private lazy val dockerConfiguration = jesConfiguration.dockerCredentials

  override def tryAbort(job: StandardAsyncJob): Unit = {
    Run(job.jobId, initializationData.genomics).abort()
  }

  override def requestsAbortAndDiesImmediately: Boolean = true

  override def receive: Receive = pollingActorClientReceive orElse super.receive

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

  private[jes] lazy val backendEngineFunctions = new JesExpressionFunctions(List(jesCallPaths.gcsPathBuilder), callContext)

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

    val writeFunctionFiles = call.task.evaluateFilesFromCommand(jobDescriptor.fullyQualifiedInputs, backendEngineFunctions) map {
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
    val globName = backendEngineFunctions.globName(wdlFile.value)
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

  override lazy val remoteStdErrPath: Path = jesStderrFile

  override lazy val remoteReturnCodePath: Path = returnCodeGcsPath

  override lazy val failOnStdErr: Boolean = runtimeAttributes.failOnStderr

  override lazy val continueOnReturnCode: ContinueOnReturnCode = runtimeAttributes.continueOnReturnCode

  override lazy val commandLineFunctions: WdlFunctions[WdlValue] = backendEngineFunctions

  override lazy val commandLinePreProcessor: (EvaluatedTaskInputs) => Try[EvaluatedTaskInputs] = mapGcsValues

  def mapGcsValues(inputs: EvaluatedTaskInputs): Try[EvaluatedTaskInputs] = {
    Try(inputs mapValues gcsPathToLocal)
  }

  override def commandLineValueMapper: (WdlValue) => WdlValue = gcsPathToLocal

  private def uploadCommandScript(command: String, withMonitoring: Boolean, globFiles: Set[WdlGlobFile]): Unit = {
    val monitoring = if (withMonitoring) {
      s"""|touch $JesMonitoringLogFile
          |chmod u+x $JesMonitoringScript
          |$JesMonitoringScript > $JesMonitoringLogFile &""".stripMargin
    } else ""

    val tmpDir = File(JesWorkingDisk.MountPoint)./("tmp").path
    val rcPath = File(JesWorkingDisk.MountPoint)./(returnCodeFilename).path
    val rcTmpPath = pathPlusSuffix(rcPath, "tmp").path

    def globManipulation(globFile: WdlGlobFile) = {

      val globDir = backendEngineFunctions.globName(globFile.value)
      val (_, disk) = relativePathAndAttachedDisk(globFile.value, runtimeAttributes.disks)
      val globDirectory = File(disk.mountPoint)./(globDir)
      val globList = File(disk.mountPoint)./(s"$globDir.list")

      s"""|mkdir $globDirectory
          |( ln -L ${globFile.value} $globDirectory 2> /dev/null ) || ( ln ${globFile.value} $globDirectory )
          |ls -1 $globDirectory > $globList
          |""".stripMargin
    }

    val globManipulations = globFiles.map(globManipulation).mkString("\n")

    val fileContent =
      s"""|#!/bin/bash
          |export _JAVA_OPTIONS=-Djava.io.tmpdir=$tmpDir
          |export TMPDIR=$tmpDir
          |$monitoring
          |(
          |cd ${JesWorkingDisk.MountPoint}
          |INSTANTIATED_COMMAND
          |)
          |echo $$? > $rcTmpPath
          |(
          |cd ${JesWorkingDisk.MountPoint}
          |$globManipulations
          |)
          |mv $rcTmpPath $rcPath
          |""".stripMargin.replace("INSTANTIATED_COMMAND", command)

    File(jesCallPaths.script).write(fileContent)
    ()
  }

  private def googleProject(descriptor: BackendWorkflowDescriptor): String = {
    descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleProject, jesAttributes.project)
  }

  private def computeServiceAccount(descriptor: BackendWorkflowDescriptor): String = {
    descriptor.workflowOptions.getOrElse(WorkflowOptionKeys.GoogleComputeServiceAccount, jesAttributes.computeServiceAccount)
  }

  override def isTerminal(runStatus: RunStatus): Boolean = {
    runStatus match {
      case _: TerminalRunStatus => true
      case _ => false
    }
  }

  private def createJesRun(jesParameters: Seq[JesParameter], runIdForResumption: Option[String]): Run = {
    Run(
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
    )
  }

  override def isFatal(throwable: Throwable): Boolean = isFatalJesException(throwable)

  override def isTransient(throwable: Throwable): Boolean = isTransientJesException(throwable)

  override def execute(): ExecutionHandle = {
    runWithJes(None)
  }

  protected def runWithJes(runIdForResumption: Option[String]): ExecutionHandle = {
    // Force runtimeAttributes to evaluate so we can fail quickly now if we need to:
    Try(runtimeAttributes) match {
      case Success(_) =>
        val command = instantiatedCommand
        val jesInputs: Set[JesInput] = generateJesInputs(jobDescriptor) ++ monitoringScript + cmdInput
        val jesOutputs: Set[JesFileOutput] = generateJesOutputs(jobDescriptor) ++ monitoringOutput
        val withMonitoring = monitoringOutput.isDefined

        val jesParameters = standardParameters ++ gcsAuthParameter ++ jesInputs ++ jesOutputs

        uploadCommandScript(command, withMonitoring, backendEngineFunctions.findGlobOutputs(call, jobDescriptor))
        val run = createJesRun(jesParameters, runIdForResumption)
        PendingExecutionHandle(jobDescriptor, StandardAsyncJob(run.runId), Option(run), previousStatus = None)
      case Failure(e) => FailedNonRetryableExecutionHandle(e)
    }
  }

  override def pollStatusAsync(handle: JesPendingExecutionHandle)
                              (implicit ec: ExecutionContext): Future[RunStatus] = {
    super[JesPollingActorClient].pollStatus(handle.runInfo.get)
  }


  override def customPollStatusFailure: PartialFunction[(ExecutionHandle, Exception), ExecutionHandle] = {
    case (oldHandle: JesPendingExecutionHandle@unchecked, e: GoogleJsonResponseException) if e.getStatusCode == 404 =>
      jobLogger.error(s"$tag JES Job ID ${oldHandle.runInfo.get.runId} has not been found, failing call")
      FailedNonRetryableExecutionHandle(e)
  }

  override lazy val startMetadataKeyValues: Map[String, Any] = super[JesJobCachingActorHelper].startMetadataKeyValues

  override def getTerminalMetadata(runStatus: RunStatus): Map[String, Any] = {
    runStatus match {
      case terminalRunStatus: TerminalRunStatus =>
        Map(
          JesMetadataKeys.MachineType -> terminalRunStatus.machineType.getOrElse("unknown"),
          JesMetadataKeys.InstanceName -> terminalRunStatus.instanceName.getOrElse("unknown"),
          JesMetadataKeys.Zone -> terminalRunStatus.zone.getOrElse("unknown")
        )
      case unknown => throw new RuntimeException(s"Attempt to get terminal metadata from non terminal status: $unknown")
    }
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
      backendEngineFunctions,
      (wdlValueToSuccess _).compose(wdlValueToGcsPath(generateJesOutputs(jobDescriptor)))
    )
  }

  override def isSuccess(runStatus: RunStatus): Boolean = {
    runStatus match {
      case _: RunStatus.Success => true
      case _: RunStatus.Failed => false
      case unknown =>
        throw new RuntimeException("isSuccess not called with RunStatus.Success or RunStatus.Failed. " +
          s"Instead got $unknown")
    }
  }

  override def handleExecutionSuccess(runStatus: RunStatus,
                                      handle: StandardAsyncPendingExecutionHandle,
                                      returnCode: Int): ExecutionHandle = {
    val success = runStatus match {
      case successStatus: RunStatus.Success => successStatus
      case unknown =>
        throw new RuntimeException(s"handleExecutionSuccess not called with RunStatus.Success. Instead got $unknown")
    }
    val outputMappings = postProcess
    val jobDetritusFiles = jesCallPaths.detritusPaths
    val executionHandle = handle
    val events = success.eventList
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

  override def handleExecutionFailure(runStatus: RunStatus,
                                      handle: StandardAsyncPendingExecutionHandle): ExecutionHandle = {
    val failed = runStatus match {
      case failedStatus: RunStatus.Failed => failedStatus
      case unknown =>
        throw new RuntimeException(s"handleExecutionFailure not called with RunStatus.Failed. Instead got $unknown")
    }
    val errorCode = failed.errorCode
    val errorMessage = failed.errorMessage

    import lenthall.numeric.IntegerUtil._

    val taskName = s"${workflowDescriptor.id}:${call.unqualifiedName}"
    val attempt = jobDescriptor.key.attempt

    if (errorMessage.exists(_.contains("Operation canceled at")))  {
      AbortedExecutionHandle
    } else if (preempted(errorCode, errorMessage)) {
      val preemptedMsg = s"Task $taskName was preempted for the ${attempt.toOrdinal} time."

      if (attempt < maxPreemption) {
        val e = PreemptedException(
          s"""$preemptedMsg The call will be restarted with another preemptible VM (max preemptible attempts number is $maxPreemption).
             |Error code $errorCode. Message: $errorMessage""".stripMargin
        )
        FailedRetryableExecutionHandle(e, None)
      } else {
        val e = PreemptedException(
          s"""$preemptedMsg The maximum number of preemptible attempts ($maxPreemption) has been reached. The call will be restarted with a non-preemptible VM.
             |Error code $errorCode. Message: $errorMessage)""".stripMargin)
        FailedRetryableExecutionHandle(e, None)
      }
    } else {
      val id = workflowDescriptor.id
      val name = jobDescriptor.call.unqualifiedName
      val message = if (errorMessage.isEmpty) "null" else errorMessage.mkString(", ")
      val exception = new RuntimeException(s"Task $id:$name failed: error code $errorCode. Message: $message")
      FailedNonRetryableExecutionHandle(exception, None)
    }
  }

  // TODO: Adapter for left over test code. Not used by main.
  private[jes] def executionResult(status: TerminalRunStatus, handle: JesPendingExecutionHandle)
                                  (implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    Future.fromTry(Try(handleExecutionResult(status, handle)))
  }

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
}
