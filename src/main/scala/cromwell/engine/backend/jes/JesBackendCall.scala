package cromwell.engine.backend.jes

import cromwell.binding._
import cromwell.binding.values.WdlValue
import cromwell.engine.AbortFunctionRegistration
import cromwell.engine.backend.BackendCall
import cromwell.util.google.GoogleCloudStoragePath

import scala.util.Try

case class JesBackendCall(backend: JesBackend,
                          workflowDescriptor: WorkflowDescriptor,
                          call: Call,
                          locallyQualifiedInputs: CallInputs,
                          callAbortRegistrationFunction: AbortFunctionRegistration) extends BackendCall {
  val callGcsPath = JesBackend.callGcsPath(workflowDescriptor.id.toString, workflowDescriptor.name, call.name)
  val callDir = GoogleCloudStoragePath(callGcsPath)
  val jesConnection = JesBackend.JesConnection
  val engineFunctions = new JesEngineFunctions(this)
  val lookupFunction: String => WdlValue = inputName => locallyQualifiedInputs.get(inputName).get

  def instantiateCommand: Try[String] = {
    val backendInputs = backend.adjustInputPaths(call, locallyQualifiedInputs)
    call.instantiateCommandLine(backendInputs, engineFunctions)
  }

  def execute: Try[Map[String, WdlValue]] = backend.execute(this)
}
