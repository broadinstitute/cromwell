package cromwell.engine.backend.jes

import java.nio.file.Paths

import com.typesafe.scalalogging.LazyLogging
import cromwell.binding._
import cromwell.binding.values.WdlFile
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.backend.jes.Run.TerminalRunStatus
import cromwell.engine.backend.jes.authentication.ProductionJesAuthentication
import cromwell.engine.backend.{BackendCall, JobKey, StdoutStderr, _}
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor}
import cromwell.util.StringDigestion._
import cromwell.util.google.GoogleCloudStoragePath

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

object JesBackendCall {
  
  def stdoutStderr(callGcsPath: String): StdoutStderr = {
    StdoutStderr(
      stdout = WdlFile(s"$callGcsPath/$StdoutFilename"),
      stderr = WdlFile(s"$callGcsPath/$StderrFilename")
    )
  }

  val StdoutFilename = "job.stdout.txt"
  val StderrFilename = "job.stderr.txt"
  val RcFilename = "job.rc.txt"

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

  val callGcsPath = backend.callGcsPath(workflowDescriptor, call.name, key.index)
  val callDir = GoogleCloudStoragePath(callGcsPath)
  val gcsExecPath = GoogleCloudStoragePath(callGcsPath + "/" + JesExecScript)

  val engineFunctions = new JesEngineFunctions(this)

  lazy val stderrJesOutput = jesOutput(callGcsPath, StderrFilename)
  lazy val stdoutJesOutput = jesOutput(callGcsPath, StdoutFilename)
  lazy val rcJesOutput = jesOutput(callGcsPath, RcFilename)
  lazy val cmdInput = JesInput(ExecParamName, gcsExecPath.toString, localFilePathFromRelativePath(JesExecScript))
  lazy val diskInput = JesInput(WorkingDiskParamName, LocalWorkingDiskValue, Paths.get(JesCromwellRoot))
  
  def standardParameters = Seq(stderrJesOutput, stdoutJesOutput, rcJesOutput, diskInput)

  def downloadRcFile = authenticated { connection => GoogleCloudStoragePath.parse(callGcsPath + "/" + RcFilename).map(connection.storage.slurpFile) }

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
          case Success(s) => handle.copy(previousStatus = Option(s)) // Copy the current handle with updated previous status.
          case Failure(e: Exception) =>
            // Log exceptions and return the original handle to try again.
            logger.warn("Caught exception, retrying: " + e.getMessage, e)
            handle
          case Failure(e: Error) => throw e // JVM-ending calamity.
          case Failure(throwable) =>
            // Someone has subclassed Throwable directly?  Beware the jackwagons, fail the execution/workflow.
            FailedExecutionHandle(throwable)
        }
      case badHandle => throw new IllegalArgumentException(s"Unexpected execution handle: $badHandle")
    }
  }

  override def resume(jobKey: JobKey)(implicit ec: ExecutionContext) = backend.resume(this, jobKey)
}
