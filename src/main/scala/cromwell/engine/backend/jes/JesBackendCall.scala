package cromwell.engine.backend.jes

import cromwell.binding._
import cromwell.engine.backend.{BackendCall, ExecutionResult}
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor}
import cromwell.util.google.GoogleCloudStoragePath

case class JesBackendCall(backend: JesBackend,
                          workflowDescriptor: WorkflowDescriptor,
                          key: CallKey,
                          locallyQualifiedInputs: CallInputs,
                          callAbortRegistrationFunction: AbortRegistrationFunction) extends BackendCall {
  val callGcsPath = backend.callGcsPath(workflowDescriptor, call.name, key.index)
  val callDir = GoogleCloudStoragePath(callGcsPath)
  val gcsExecPath = GoogleCloudStoragePath(callGcsPath + "/exec.sh")
  val jesConnection = backend.JesConnection
  val engineFunctions = new JesEngineFunctions(this)
  def execute: ExecutionResult = backend.execute(this)
}
