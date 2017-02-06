package cromwell.backend.impl.jes

import java.net.SocketTimeoutException
import java.nio.file.{Path, Paths}

import akka.actor.ActorRef
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import cromwell.backend._
import cromwell.backend.async.{AbortedExecutionHandle, ExecutionHandle, FailedNonRetryableExecutionHandle, FailedRetryableExecutionHandle, PendingExecutionHandle}
import cromwell.backend.impl.jes.RunStatus.TerminalRunStatus
import cromwell.backend.impl.jes.io._
import cromwell.backend.impl.jes.statuspolling.JesPollingActorClient
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams, StandardAsyncJob}
import cromwell.core._
import cromwell.core.path.PathImplicits._
import cromwell.core.path.proxy.PathProxy
import cromwell.core.retry.SimpleExponentialBackoff
import wdl4s._
import wdl4s.expression.NoFunctions
import wdl4s.values._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object JesAsyncBackendJobExecutionActor {
  val JesOperationIdKey = "__jes_operation_id"

  object WorkflowOptionKeys {
    val MonitoringScript = "monitoring_script"
    val GoogleProject = "google_project"
    val GoogleComputeServiceAccount = "google_compute_service_account"
  }

  type JesPendingExecutionHandle =
    PendingExecutionHandle[StandardAsyncJob, Run, RunStatus]

  private val ExtraConfigParamName = "__extra_config_gcs_path"
}

class JesAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends BackendJobLifecycleActor with StandardAsyncExecutionActor with JesJobCachingActorHelper
    with JesPollingActorClient {

  import JesAsyncBackendJobExecutionActor._

  val jesBackendSingletonActor: ActorRef =
    standardParams.backendSingletonActorOption.getOrElse(
      throw new RuntimeException("JES Backend actor cannot exist without the JES backend singleton actor"))

  override type StandardAsyncRunInfo = Run

  override type StandardAsyncRunStatus = RunStatus

  override val pollingActor: ActorRef = jesBackendSingletonActor

  override lazy val pollBackOff = SimpleExponentialBackoff(
    initialInterval = 30 seconds, maxInterval = jesAttributes.maxPollingInterval seconds, multiplier = 1.1)

  override lazy val executeOrRecoverBackOff = SimpleExponentialBackoff(
    initialInterval = 3 seconds, maxInterval = 20 seconds, multiplier = 1.1)

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

  override lazy val commandDirectory: Path = JesWorkingDisk.MountPoint

  override def commandScriptPreamble: String = {
    if (monitoringOutput.isDefined) {
      s"""|touch $JesMonitoringLogFile
          |chmod u+x $JesMonitoringScript
          |$JesMonitoringScript > $JesMonitoringLogFile &""".stripMargin
    } else ""
  }

  override def globParentDirectory(wdlGlobFile: WdlGlobFile): Path = {
    val (_, disk) = relativePathAndAttachedDisk(wdlGlobFile.value, runtimeAttributes.disks)
    disk.mountPoint
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
        val jesInputs: Set[JesInput] = generateJesInputs(jobDescriptor) ++ monitoringScript + cmdInput
        val jesOutputs: Set[JesFileOutput] = generateJesOutputs(jobDescriptor) ++ monitoringOutput

        val jesParameters = standardParameters ++ gcsAuthParameter ++ jesInputs ++ jesOutputs

        jobPaths.script.writeAsText(commandScriptContents)
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

  override def mapOutputWdlFile(wdlFile: WdlFile): WdlFile = {
    wdlFileToGcsPath(generateJesOutputs(jobDescriptor))(wdlFile)
  }

  private[jes] def wdlFileToGcsPath(jesOutputs: Set[JesFileOutput])(wdlFile: WdlFile): WdlFile = {
    jesOutputs collectFirst {
      case jesOutput if jesOutput.name == makeSafeJesReferenceName(wdlFile.valueString) => WdlFile(jesOutput.gcs)
    } getOrElse wdlFile
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

  override def getTerminalEvents(runStatus: RunStatus): Seq[ExecutionEvent] = {
    runStatus match {
      case successStatus: RunStatus.Success => successStatus.eventList
      case unknown =>
        throw new RuntimeException(s"handleExecutionSuccess not called with RunStatus.Success. Instead got $unknown")
    }
  }

  override def retryEvaluateOutputs(exception: Exception): Boolean = {
    exception match {
      case aggregated: CromwellAggregatedException =>
        aggregated.throwables.collectFirst { case s: SocketTimeoutException => s }.isDefined
      case _ => false
    }
  }

  private def extractErrorCodeFromErrorMessage(errorMessage: String): Int = {
    errorMessage.substring(0, errorMessage.indexOf(':')).toInt
  }

  private def preempted(errorCode: Int, errorMessage: String): Boolean = {
    def isPreemptionCode(code: Int) = code == 13 || code == 14

    try {
      errorCode == 10 && isPreemptionCode(extractErrorCodeFromErrorMessage(errorMessage)) && preemptible
    } catch {
      case _: NumberFormatException | _: StringIndexOutOfBoundsException =>
        jobLogger.warn(s"Unable to parse JES error code from error messages: [{}], assuming this was not a preempted VM.", errorMessage.mkString(", "))
        false
    }
  }

  override def handleExecutionFailure(runStatus: RunStatus,
                                      handle: StandardAsyncPendingExecutionHandle,
                                      returnCode: Option[Int]): ExecutionHandle = {
    import lenthall.numeric.IntegerUtil._

    val failed = runStatus match {
      case failedStatus: RunStatus.Failed => failedStatus
      case unknown => throw new RuntimeException(s"handleExecutionFailure not called with RunStatus.Failed. Instead got $unknown")
    }

    // In the event that the error isn't caught by the special cases, this will be used
    def nonRetryableException = {
      val exception = failed.toFailure(jobPaths.jobKey.tag, Option(jobPaths.stderr))
      FailedNonRetryableExecutionHandle(exception, returnCode)
    }

    val errorCode = failed.errorCode
    val taskName = s"${workflowDescriptor.id}:${call.unqualifiedName}"
    val attempt = jobDescriptor.key.attempt

    failed.errorMessage match {
      case Some(errorMessage) =>
        if (errorMessage.contains("Operation canceled at"))  {
          AbortedExecutionHandle
        } else if (preempted(errorCode, errorMessage)) {
          val preemptedMsg = s"Task $taskName was preempted for the ${attempt.toOrdinal} time."
          if (attempt < maxPreemption) {
            val e = PreemptedException(
              s"""$preemptedMsg The call will be restarted with another preemptible VM (max preemptible attempts number is $maxPreemption).
                 |Error code $errorCode. Message: $errorMessage""".stripMargin
            )
            FailedRetryableExecutionHandle(e, returnCode)
          } else {
            val e = PreemptedException(
              s"""$preemptedMsg The maximum number of preemptible attempts ($maxPreemption) has been reached. The call will be restarted with a non-preemptible VM.
                 |Error code $errorCode. Message: $errorMessage)""".stripMargin)
            FailedRetryableExecutionHandle(e, returnCode)
          }
        } else { nonRetryableException }
      case None => nonRetryableException
    }
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

  override def mapCommandLineWdlFile(wdlFile: WdlFile): WdlFile = {
    getPath(wdlFile.valueString) match {
      case Success(gcsPath: CloudStoragePath) =>
        WdlFile(localFilePathFromCloudStoragePath(workingDisk.mountPoint, gcsPath).toString, wdlFile.isGlob)
      case Success(proxy: PathProxy) =>
        proxy.unbox(classOf[CloudStoragePath]) map { gcsPath =>
          WdlFile(localFilePathFromCloudStoragePath(workingDisk.mountPoint, gcsPath).toString, wdlFile.isGlob)
        } getOrElse wdlFile
      case _ => wdlFile
    }
  }
}
