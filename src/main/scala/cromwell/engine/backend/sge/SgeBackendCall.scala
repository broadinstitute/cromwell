package cromwell.engine.backend.sge

import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, CallInputs, WorkflowDescriptor}
import cromwell.engine.AbortFunctionRegistration
import cromwell.engine.backend.{LocalFileSystemBackendCall, BackendCall}
import cromwell.engine.backend.local.LocalBackend

import scala.util.Try

case class SgeBackendCall(backend: SgeBackend,
                          workflowDescriptor: WorkflowDescriptor,
                          call: Call,
                          locallyQualifiedInputs: CallInputs,
                          callAbortRegistrationFunction: AbortFunctionRegistration) extends BackendCall with LocalFileSystemBackendCall {
  val workflowRootPath = LocalBackend.hostExecutionPath(workflowDescriptor)
  val callRootPath = LocalBackend.hostCallPath(workflowDescriptor, call.name)
  val stdout = callRootPath.resolve("stdout")
  val stderr = callRootPath.resolve("stderr")
  val script = callRootPath.resolve("script.sh")
  val rc = callRootPath.resolve("rc")
  val engineFunctions: SgeEngineFunctions = new SgeEngineFunctions(callRootPath, stdout, stderr)
  val lookupFunction: String => WdlValue = inputName => locallyQualifiedInputs.get(inputName).get
  callRootPath.toFile.mkdirs

  def instantiateCommand: Try[String] = {
    val backendInputs = backend.adjustInputPaths(call, locallyQualifiedInputs)
    call.instantiateCommandLine(backendInputs, engineFunctions)
  }

  def execute: Try[Map[String, WdlValue]] = backend.execute(this)
}
