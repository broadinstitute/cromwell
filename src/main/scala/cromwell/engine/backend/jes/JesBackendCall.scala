package cromwell.engine.backend.jes

import java.io.File
import java.nio.file.Paths

import cromwell.binding._
import cromwell.binding.values.WdlFile
import cromwell.engine.backend.jes.JesBackend.{JesInput, JesOutput}
import cromwell.engine.backend.jes.authentication.ProductionJesAuthentication
import cromwell.engine.backend.{BackendCall, ExecutionResult, StdoutStderr}
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor}
import cromwell.util.google.GoogleCloudStoragePath


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

class JesBackendCall(val backend: JesBackend,
                     val workflowDescriptor: WorkflowDescriptor,
                     val key: CallKey,
                     val locallyQualifiedInputs: CallInputs,
                     val callAbortRegistrationFunction: AbortRegistrationFunction) extends BackendCall with ProductionJesAuthentication {

  import JesBackendCall._

  def jesCommandLine = s"/bin/bash exec.sh > $StdoutFilename 2> $StderrFilename"

  val callGcsPath = backend.callGcsPath(workflowDescriptor, call.name, key.index)
  val callDir = GoogleCloudStoragePath(callGcsPath)
  val gcsExecPath = GoogleCloudStoragePath(callGcsPath + "/exec.sh")

  val engineFunctions = new JesEngineFunctions(this)
  
  lazy val stderrJesOutput = jesOutput(callGcsPath, StderrFilename)
  lazy val stdoutJesOutput = jesOutput(callGcsPath, StdoutFilename)
  lazy val rcJesOutput = jesOutput(callGcsPath, RcFilename)
  lazy val diskInput = JesInput("working_disk", "disk://local-disk", new File(JesBackend.JesCromwellRoot).toPath)
  
  def standardParameters = Seq(stderrJesOutput, stdoutJesOutput, rcJesOutput, diskInput)
  def rcGcsPath = rcJesOutput.gcs

  def execute: ExecutionResult = backend.execute(this)

  def downloadRcFile = authenticated { connection => GoogleCloudStoragePath.parse(callGcsPath + "/" + RcFilename).map(connection.storage.slurpFile) }

}
