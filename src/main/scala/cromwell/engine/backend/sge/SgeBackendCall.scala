package cromwell.engine.backend.sge

import cromwell.binding.CallInputs
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.backend._
import cromwell.engine.backend.{JobKey, BackendCall, ExecutionResult, LocalFileSystemBackendCall}
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor}

import scala.concurrent.Future

case class SgeBackendCall(backend: SgeBackend,
                          workflowDescriptor: WorkflowDescriptor,
                          key: CallKey,
                          locallyQualifiedInputs: CallInputs,
                          callAbortRegistrationFunction: AbortRegistrationFunction) extends BackendCall with LocalFileSystemBackendCall {
  val workflowRootPath = LocalBackend.hostExecutionPath(workflowDescriptor)
  val callRootPath = LocalBackend.hostCallPath(workflowDescriptor, call.name, key.index)
  val stdout = callRootPath.resolve("stdout")
  val stderr = callRootPath.resolve("stderr")
  val script = callRootPath.resolve("script.sh")
  val returnCode = callRootPath.resolve("rc")
  val engineFunctions: SgeEngineFunctions = new SgeEngineFunctions(callRootPath, stdout, stderr)
  callRootPath.toFile.mkdirs

  // Block for the execution.  TODO don't block.
  override def execute = Future.successful(CompletedExecutionHandle(backend.execute(this)))

  override def resume(jobKey: JobKey) = Future.failed(new NotImplementedError)

  override def poll(previous: ExecutionHandle) = Future.successful(previous)
}
