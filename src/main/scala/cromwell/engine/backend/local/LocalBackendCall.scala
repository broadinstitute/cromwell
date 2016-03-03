package cromwell.engine.backend.local

import better.files._
import cromwell.engine.backend.{BackendCall, LocalFileSystemBackendCall, _}
import cromwell.engine.workflow.BackendCallKey
import cromwell.engine.{AbortRegistrationFunction, CallContext, WorkflowDescriptor}
import wdl4s.CallInputs

import scala.concurrent.{ExecutionContext, Future}

case class LocalBackendCall(backend: LocalBackend,
                            workflowDescriptor: WorkflowDescriptor,
                            key: BackendCallKey,
                            locallyQualifiedInputs: CallInputs,
                            callAbortRegistrationFunction: Option[AbortRegistrationFunction]) extends BackendCall with LocalFileSystemBackendCall {
  val dockerContainerExecutionDir = workflowDescriptor.workflowRootPathWithBaseRoot(LocalBackend.ContainerRoot)
  lazy val containerCallRoot = runtimeAttributes.docker match {
    case Some(docker) => callRootPathWithBaseRoot(LocalBackend.ContainerRoot)
    case None => callRootPath
  }
  val returnCode = callRootPath.resolve("rc")
  val stdout = callRootPath.resolve("stdout")
  val stderr = callRootPath.resolve("stderr")
  val script = callRootPath.resolve("script")
  private val callContext = new CallContext(callRootPath.fullPath, stdout.fullPath, stderr.fullPath)
  val engineFunctions = new LocalCallEngineFunctions(workflowDescriptor.ioManager, callContext)

  callRootPath.toFile.mkdirs

  override def execute(implicit ec: ExecutionContext) = backend.execute(this)

  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext) = Future.successful(previous)

  override def useCachedCall(avoidedTo: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] =
    backend.useCachedCall(avoidedTo.asInstanceOf[LocalBackendCall], this)

  override def stdoutStderr: CallLogs = backend.stdoutStderr(this)
}
