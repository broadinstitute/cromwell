package cromwell.engine.backend.local

import cromwell.binding.values.WdlValue
import cromwell.binding.{CallInputs, WorkflowDescriptor}
import cromwell.engine.AbortRegistrationFunction
import cromwell.engine.backend.{BackendCall, ExecutionResult, LocalFileSystemBackendCall}
import cromwell.engine.workflow.CallKey

case class LocalBackendCall(backend: LocalBackend,
                            workflowDescriptor: WorkflowDescriptor,
                            key: CallKey,
                            locallyQualifiedInputs: CallInputs,
                            callAbortRegistrationFunction: AbortRegistrationFunction) extends BackendCall with LocalFileSystemBackendCall {
  val workflowRootPath = LocalBackend.hostExecutionPath(workflowDescriptor)
  val callRootPath = LocalBackend.hostCallPath(workflowDescriptor, call.name, key.index)
  val dockerContainerExecutionDir = LocalBackend.containerExecutionPath(workflowDescriptor)
  val containerCallRoot = call.docker match {
    case Some(docker) => LocalBackend.containerCallPath(workflowDescriptor, call.name, key.index)
    case None => callRootPath
  }
  val rc = callRootPath.resolve("rc")
  val stdout = callRootPath.resolve("stdout")
  val stderr = callRootPath.resolve("stderr")
  val script = callRootPath.resolve("script")
  val engineFunctions: LocalEngineFunctions = new LocalEngineFunctions(callRootPath, stdout, stderr)
  val lookupFunction: String => WdlValue = inputName => locallyQualifiedInputs.get(inputName).get
  callRootPath.toFile.mkdirs
  def execute: ExecutionResult = backend.execute(this)
}
