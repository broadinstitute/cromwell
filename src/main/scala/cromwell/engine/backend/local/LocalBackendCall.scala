package cromwell.engine.backend.local

import cromwell.binding.CallInputs
import cromwell.engine.backend._
import cromwell.engine.backend.{JobKey, BackendCall, ExecutionResult, LocalFileSystemBackendCall}
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor}

import scala.concurrent.Future

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
  val returnCode = callRootPath.resolve("rc")
  val stdout = callRootPath.resolve("stdout")
  val stderr = callRootPath.resolve("stderr")
  val script = callRootPath.resolve("script")
  val engineFunctions: LocalEngineFunctions = new LocalEngineFunctions(callRootPath, stdout, stderr)
  callRootPath.toFile.mkdirs

  // Block for the execution.  TODO don't block.
  override def execute = Future.successful(CompletedExecutionHandle(backend.execute(this)))

  override def resume(jobKey: JobKey) = Future.failed(new NotImplementedError)

  override def poll(previous: ExecutionHandle) = Future.successful(previous)
}
