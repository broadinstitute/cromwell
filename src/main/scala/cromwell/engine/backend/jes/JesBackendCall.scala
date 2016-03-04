package cromwell.engine.backend.jes

import java.nio.file.Paths

import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.typesafe.scalalogging.LazyLogging
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.backend.jes.Run.TerminalRunStatus
import cromwell.engine.backend.jes.authentication.ProductionJesAuthentication
import cromwell.engine.backend.{BackendCall, CallLogs, JobKey, _}
import cromwell.engine.io.gcs.GcsPath
import cromwell.engine.workflow.BackendCallKey
import cromwell.engine.{AbortRegistrationFunction, CallContext, WorkflowDescriptor}
import wdl4s._
import wdl4s.values._

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object JesBackendCall {
  def jesLogBasename(key: BackendCallKey) = {
    val index = key.index.map(s => s"-$s").getOrElse("")
    s"${key.scope.unqualifiedName}$index"
  }

  def jesLogFilename(key: BackendCallKey) = s"${jesLogBasename(key)}.log"
  def jesLogStdoutFilename(key: BackendCallKey) = s"${jesLogBasename(key)}-stdout.log"
  def jesLogStderrFilename(key: BackendCallKey) = s"${jesLogBasename(key)}-stderr.log"
  def jesReturnCodeFilename(key: BackendCallKey) = s"${jesLogBasename(key)}-rc.txt"
}

class JesBackendCall(val backend: JesBackend,
                     val jobDescriptor: BackendCallJobDescriptor,
                     val callAbortRegistrationFunction: Option[AbortRegistrationFunction])
  extends BackendCall with ProductionJesAuthentication with LazyLogging {

  import JesBackend._
  import JesBackendCall._

  // TODO: Assuming that runtimeAttributes.disks always has a 'local-disk'
  lazy val workingDisk = runtimeAttributes.disks.find(_.name == JesWorkingDisk.Name).get

  def jesCommandLine = s"/bin/bash ${cmdInput.containerPath.toAbsolutePath.toString}"

  lazy val callGcsPath = key.callRootPathWithBaseRoot(workflowDescriptor, backend.rootPath(workflowDescriptor.workflowOptions))
  //TODO: Switch to NioGcsPath
  lazy val callDir = GcsPath(callGcsPath.toString)
  lazy val gcsExecPath = GcsPath(callGcsPath + "/" + JesExecScript)
  lazy val defaultMonitoringOutputPath = callGcsPath + "/" + JesMonitoringLogFile

  lazy val returnCodeFilename = jesReturnCodeFilename(key)
  lazy val jesStdoutGcsPath = s"$callGcsPath/${jesLogStdoutFilename(key)}"
  lazy val jesStderrGcsPath = s"$callGcsPath/${jesLogStderrFilename(key)}"
  lazy val jesLogGcsPath = s"$callGcsPath/${jesLogFilename(key)}"
  lazy val returnCodeGcsPath = s"$callGcsPath/$returnCodeFilename"

  private lazy val callContext = new CallContext(callGcsPath.toString, jesStdoutGcsPath, jesStderrGcsPath)

  lazy val callEngineFunctions = new JesCallEngineFunctions(workflowDescriptor.ioManager, callContext)

  lazy val rcJesOutput = JesFileOutput(returnCodeFilename, returnCodeGcsPath, Paths.get(returnCodeFilename), workingDisk)
  lazy val cmdInput = JesFileInput(ExecParamName, gcsExecPath.toString, Paths.get(JesExecScript), workingDisk)

  def standardParameters = Seq(rcJesOutput)

  def downloadRcFile = authenticateAsUser(workflowDescriptor) { storage => Try(storage.readFile(returnCodeGcsPath)) }

  def instantiateCommand: Try[String] = {
    val backendInputs = backend.adjustInputPaths(this)
    call.instantiateCommandLine(backendInputs, callEngineFunctions, JesBackend.gcsPathToLocal)
  }

  /**
    * Determines the maximum number of times a call can be started with a Preemptible VM.
    * TODO: Use configuration as a way to set this globally.
    * Currently workflow options act as default for runtime attributes, configuration could do the same for workflow options.
    */
  lazy val maxPreemption = runtimeAttributes.preemptible

  lazy val preemptible = key.attempt <= maxPreemption

  /**
   * Determine the output directory for the files matching a particular glob.
   */
  def globOutputPath(glob: String) = s"$callGcsPath/glob-${glob.md5Sum}/"

  override def execute(implicit ec: ExecutionContext) = backend.execute(this)

  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext) = Future {
    previous match {
      case handle: JesPendingExecutionHandle =>
        val status = Try(handle.run.checkStatus(this, handle.previousStatus))
        status match {
          case Success(s: TerminalRunStatus) => backend.executionResult(s, handle)
          case Success(s) => handle.copy(previousStatus = Option(s)).future // Copy the current handle with updated previous status.
          case Failure(e: GoogleJsonResponseException) if e.getStatusCode == 404 =>
            logger.error(s"JES Job ID ${handle.run.runId} has not been found, failing call")
            FailedExecutionHandle(e).future
          case Failure(e: Exception) =>
            // Log exceptions and return the original handle to try again.
            logger.warn("Caught exception, retrying: " + e.getMessage, e)
            handle.future
          case Failure(e: Error) => Future.failed(e) // JVM-ending calamity.
          case Failure(throwable) =>
            // Someone has subclassed Throwable directly?
            FailedExecutionHandle(throwable).future
        }
      case f: FailedExecutionHandle => f.future
      case s: SuccessfulExecutionHandle => s.future
      case badHandle => Future.failed(new IllegalArgumentException(s"Unexpected execution handle: $badHandle"))
    }
  } flatten

  override def resume(jobKey: JobKey)(implicit ec: ExecutionContext) = backend.resume(this, jobKey)

  override def useCachedCall(avoidedTo: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] =
    backend.useCachedCall(avoidedTo.asInstanceOf[JesBackendCall], this)

  override def stdoutStderr: CallLogs = {
    CallLogs(
      stdout = WdlFile(jesStdoutGcsPath),
      stderr = WdlFile(jesStderrGcsPath),
      Option(Map("log" -> WdlFile(jesLogGcsPath)))
    )
  }
}
