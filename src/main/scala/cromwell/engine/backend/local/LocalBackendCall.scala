package cromwell.engine.backend.local

import cromwell.binding.values.WdlValue
import cromwell.binding.{CallOutputs, CallInputs, Call, WorkflowDescriptor}
import cromwell.engine.AbortFunctionRegistration
import cromwell.engine.backend.{LocalFileSystemBackendCall, BackendCall}

import scala.util.Try

case class LocalBackendCall(backend: LocalBackend,
                            workflowDescriptor: WorkflowDescriptor,
                            call: Call,
                            locallyQualifiedInputs: CallInputs,
                            callAbortRegistrationFunction: AbortFunctionRegistration) extends BackendCall with LocalFileSystemBackendCall {
  val workflowRootPath = LocalBackend.hostExecutionPath(workflowDescriptor)
  val callRootPath = LocalBackend.hostCallPath(workflowDescriptor, call.name)
  val stdout = callRootPath.resolve("stdout")
  val stderr = callRootPath.resolve("stderr")
  val script = callRootPath.resolve("script")
  val engineFunctions: LocalEngineFunctions = new LocalEngineFunctions(callRootPath, stdout, stderr)
  val lookupFunction: String => WdlValue = inputName => locallyQualifiedInputs.get(inputName).get
  callRootPath.toFile.mkdirs

  def instantiateCommand: Try[String] = {
    val backendInputs = backend.adjustInputPaths(call, locallyQualifiedInputs)
    call.instantiateCommandLine(backendInputs, engineFunctions)
  }
  
  def execute: Try[CallOutputs] = backend.execute(this)
}
