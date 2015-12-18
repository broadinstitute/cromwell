package cromwell.engine.backend.sge

import cromwell.binding.CallInputs
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.backend.{BackendCall, LocalFileSystemBackendCall, _}
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor}

import scala.concurrent.{ExecutionContext, Future}

case class SgeBackendCall(backend: SgeBackend,
                          workflowDescriptor: WorkflowDescriptor,
                          key: CallKey,
                          locallyQualifiedInputs: CallInputs,
                          callAbortRegistrationFunction: AbortRegistrationFunction) extends BackendCall with LocalFileSystemBackendCall {
  val workflowRootPath = LocalBackend.hostExecutionPath(workflowDescriptor)
  val callRootPath = LocalBackend.hostCallPath(workflowDescriptor, call.unqualifiedName, key.index)
  val stdout = callRootPath.resolve("stdout")
  val stderr = callRootPath.resolve("stderr")
  val script = callRootPath.resolve("script.sh")
  val returnCode = callRootPath.resolve("rc")
  val engineFunctions: SgeEngineFunctions = new SgeEngineFunctions(callRootPath, stdout, stderr, workflowDescriptor.IOInterface)
  callRootPath.toFile.mkdirs

  override def execute(implicit ec: ExecutionContext) = backend.execute(this)

  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext) = Future.successful(previous)

  override def useCachedCall(avoidedTo: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] =
    backend.useCachedCall(avoidedTo.asInstanceOf[SgeBackendCall], this)
}
