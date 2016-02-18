package cromwell.engine.backend.sge

import cromwell.engine.backend.{BackendCall, LocalFileSystemBackendCall, _}
import cromwell.engine.workflow.BackendCallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor}
import wdl4s.CallInputs

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

case class SgeBackendCall(backend: SgeBackend,
                          workflowDescriptor: WorkflowDescriptor,
                          key: BackendCallKey,
                          locallyQualifiedInputs: CallInputs,
                          callAbortRegistrationFunction: Option[AbortRegistrationFunction]) extends BackendCall with LocalFileSystemBackendCall {
  val workflowRootPath = workflowDescriptor.workflowRootPath
  val stdout = callRootPath.resolve("stdout")
  val stderr = callRootPath.resolve("stderr")
  val script = callRootPath.resolve("script.sh")
  val returnCode = callRootPath.resolve("rc")
  val engineFunctions: SgeEngineFunctions = new SgeEngineFunctions(callRootPath, stdout, stderr, workflowDescriptor.ioManager)
  callRootPath.toFile.mkdirs

  def instantiateCommand: Try[String] = {
    val backendInputs = backend.adjustInputPaths(this)
    call.instantiateCommandLine(backendInputs, engineFunctions)
  }

  override def execute(implicit ec: ExecutionContext) = backend.execute(this)

  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext) = Future.successful(previous)

  override def useCachedCall(avoidedTo: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] =
    backend.useCachedCall(avoidedTo.asInstanceOf[SgeBackendCall], this)

  override def stdoutStderr: CallLogs = backend.stdoutStderr(this)
}
