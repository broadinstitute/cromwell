package cromwell.engine.backend.jes

import java.nio.file.Paths

import cromwell.binding._
import cromwell.binding.values.WdlFile
import cromwell.engine.backend.jes.JesBackend._
import cromwell.engine.backend.jes.Run.TerminalRunStatus
import cromwell.engine.backend.jes.authentication.ProductionJesAuthentication
import cromwell.engine.backend._
import cromwell.engine.backend.{JobKey, BackendCall, ExecutionResult, StdoutStderr}
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor}
import cromwell.util.StringDigestion._
import cromwell.util.google.GoogleCloudStoragePath

import scala.concurrent.Future

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
                     val callAbortRegistrationFunction: AbortRegistrationFunction) extends BackendCall with ProductionJesAuthentication {

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

  /** TODO creating the JES run should be async but currently is not. */
  override def execute = Future.successful(backend.execute(this))

  /** Note this is currently not async querying JES status, only in not blocking *between*
    *  JES status queries. */
  override def poll(previous: ExecutionHandle) = previous match {
    case handle: JesPendingExecutionHandle =>
      val status = handle.run.checkStatus(this, handle.previousStatus)
      val nextHandle = status match {
        case s: TerminalRunStatus => CompletedExecutionHandle(backend.executionResult(s, handle))
        case s => handle.copy(previousStatus = Option(s)) // Copy the current handle with updated previous status.
      }
      Future.successful(nextHandle)
    case badHandle => Future.failed(new IllegalArgumentException(s"Unexpected execution handle: $badHandle"))
  }

  /**
   * TODO resuming a JES run should be async but currently is not, though this wouldn't be nearly as big a win as
   * for the execute case since this probably doesn't do much of anything with JES.
   */
  override def resume(jobKey: JobKey) = Future.successful(backend.resume(this, jobKey))
}
