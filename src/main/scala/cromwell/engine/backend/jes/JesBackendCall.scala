package cromwell.engine.backend.jes

import java.nio.file.Paths

import cromwell.binding._
import cromwell.binding.values.WdlFile
import cromwell.engine.backend.jes.JesBackend.JesOutput
import cromwell.engine.backend.{StdoutStderr, BackendCall, ExecutionResult}
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor}
import cromwell.util.google.GoogleCloudStoragePath

import scala.util.Try


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
  private def jesOutput(callGcsPath: String, filename: String): JesOutput = JesOutput(filename, s"$callGcsPath/$filename", Paths.get(filename))
}

case class JesBackendCall(backend: JesBackend,
                          workflowDescriptor: WorkflowDescriptor,
                          key: CallKey,
                          locallyQualifiedInputs: CallInputs,
                          callAbortRegistrationFunction: AbortRegistrationFunction) extends BackendCall {

  import JesBackendCall._

  def jesCommandLine = s"/bin/bash exec.sh > $StdoutFilename 2> $StderrFilename"

  val callGcsPath = backend.callGcsPath(workflowDescriptor, call.name, key.index)
  val callDir = GoogleCloudStoragePath(callGcsPath)
  val gcsExecPath = GoogleCloudStoragePath(callGcsPath + "/exec.sh")
  val jesConnection = backend.JesConnection
  val engineFunctions = new JesEngineFunctions(this)
  
  lazy val stderrJesOutput = jesOutput(callGcsPath, StderrFilename)
  lazy val stdoutJesOutput = jesOutput(callGcsPath, StdoutFilename)
  lazy val rcJesOutput = jesOutput(callGcsPath, RcFilename)
  
  def standardParameters = Seq(stderrJesOutput, stdoutJesOutput, rcJesOutput) 
  def rcGcsPath = rcJesOutput.gcs
  def execute: ExecutionResult = backend.execute(this)
  def downloadRcFile: Try[String] = GoogleCloudStoragePath.parse(callGcsPath + "/" + RcFilename).map(jesConnection.storage.slurpFile)
}
