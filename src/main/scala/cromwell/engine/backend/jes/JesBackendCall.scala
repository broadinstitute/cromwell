package cromwell.engine.backend.jes

import java.nio.file.Paths

import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.typesafe.scalalogging.LazyLogging
import wdl4s._
import wdl4s.values.WdlFile
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.backend.jes.Run.TerminalRunStatus
import cromwell.engine.backend.jes.authentication.ProductionJesAuthentication
import cromwell.engine.backend.{BackendCall, CallLogs, JobKey, _}
import cromwell.engine.io.gcs.GcsPath
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, CallContext, WorkflowDescriptor, _}

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import Hashing._

object JesBackendCall {

  def stdoutStderr(callGcsPath: String): CallLogs = {
    CallLogs(
      stdout = WdlFile(s"$callGcsPath/$StdoutFilename"),
      stderr = WdlFile(s"$callGcsPath/$StderrFilename"),
      Option(Map("log" -> WdlFile(s"$callGcsPath/$JesLog"),
        "stdout" -> WdlFile(s"$callGcsPath/$JesStdout"),
        "stderr" -> WdlFile(s"$callGcsPath/$JesStderr")
      ))
    )
  }

  val StdoutFilename = "job.stdout.txt"
  val StderrFilename = "job.stderr.txt"
  val RcFilename = "job.rc.txt"

  val JesLogBasename = "jes"
  private [jes] val JesLog = s"$JesLogBasename.log"
  private val JesStdout = s"$JesLogBasename-stdout.log"
  private val JesStderr = s"$JesLogBasename-stderr.log"

  private def jesOutput(callGcsPath: String, filename: String): JesOutput =
    JesOutput(filename, s"$callGcsPath/$filename", localFilePathFromRelativePath(filename))
}

class JesBackendCall(val backend: JesBackend,
                     val workflowDescriptor: WorkflowDescriptor,
                     val key: CallKey,
                     val locallyQualifiedInputs: CallInputs,
                     val callAbortRegistrationFunction: AbortRegistrationFunction)
  extends BackendCall with ProductionJesAuthentication with LazyLogging {

  import JesBackend._
  import JesBackendCall._

  def jesCommandLine = s"/bin/bash ${cmdInput.local} > ${stdoutJesOutput.local} 2> ${stderrJesOutput.local}"

  val callGcsPath = backend.callGcsPath(workflowDescriptor, call.unqualifiedName, key.index)
  val callDir = GcsPath(callGcsPath)
  val gcsExecPath = GcsPath(callGcsPath + "/" + JesExecScript)
  val defaultMonitoringOutputPath = callGcsPath + "/" + JesMonitoringLogFile

  private val callContext = new CallContext(callGcsPath, stdoutJesOutput.gcs, stderrJesOutput.gcs)

  val engineFunctions = new JesCallEngineFunctions(workflowDescriptor.ioManager, callContext)

  lazy val stderrJesOutput = jesOutput(callGcsPath, StderrFilename)
  lazy val stdoutJesOutput = jesOutput(callGcsPath, StdoutFilename)
  lazy val rcJesOutput = jesOutput(callGcsPath, RcFilename)
  lazy val cmdInput = JesInput(ExecParamName, gcsExecPath.toString, localFilePathFromRelativePath(JesExecScript))
  lazy val diskInput = JesInput(WorkingDiskParamName, LocalWorkingDiskValue, Paths.get(JesCromwellRoot))

  def standardParameters = Seq(stderrJesOutput, stdoutJesOutput, rcJesOutput, diskInput)

  def downloadRcFile = authenticateAsUser(workflowDescriptor) { storage => Try(storage.readFile(s"$callGcsPath/$RcFilename")) }

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
}
