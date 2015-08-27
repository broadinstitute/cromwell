package cromwell.engine.backend.jes

import cromwell.binding._
import cromwell.binding.values.WdlValue
import cromwell.engine.AbortRegistrationFunction
import cromwell.engine.backend.BackendCall
import cromwell.util.google.GoogleCloudStoragePath

import scala.util.Try

case class JesBackendCall(backend: JesBackend,
                          workflowDescriptor: WorkflowDescriptor,
                          call: Call,
                          locallyQualifiedInputs: CallInputs,
                          callAbortRegistrationFunction: AbortRegistrationFunction) extends BackendCall {
  val callGcsPath = JesBackend.callGcsPath(workflowDescriptor.id.toString, workflowDescriptor.name, call.name)
  val callDir = GoogleCloudStoragePath(callGcsPath)
  val gcsExecPath = GoogleCloudStoragePath(callGcsPath + "/exec.sh")
  val jesConnection = JesBackend.JesConnection
  val engineFunctions = new JesEngineFunctions(this)
  val lookupFunction: String => WdlValue = inputName => locallyQualifiedInputs.get(inputName).get
  def execute: Try[Map[String, WdlValue]] = backend.execute(this)
}
