package cromwell.engine.backend.sge

import cromwell.engine.AbortRegistrationFunction
import cromwell.engine.backend.{BackendCall, LocalFileSystemBackendCall, _}

import scala.concurrent.{ExecutionContext, Future}

case class SgeBackendCall(backend: SgeBackend,
                          jobDescriptor: BackendCallJobDescriptor,
                          callAbortRegistrationFunction: Option[AbortRegistrationFunction]) extends BackendCall with LocalFileSystemBackendCall {

  val returnCode = backend.returnCode(jobDescriptor)

  override def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext) = Future.successful(previous)

  override def stdoutStderr: CallLogs = backend.stdoutStderr(jobDescriptor)
}
