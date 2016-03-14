package cromwell.engine.backend.jes

import com.typesafe.scalalogging.LazyLogging
import cromwell.engine.AbortRegistrationFunction
import cromwell.engine.backend.jes.authentication.ProductionJesAuthentication
import cromwell.engine.backend.{BackendCall, CallLogs, _}

import scala.language.postfixOps

class JesBackendCall(val backend: JesBackend,
                     val jobDescriptor: BackendCallJobDescriptor,
                     val callAbortRegistrationFunction: Option[AbortRegistrationFunction])
  extends BackendCall with ProductionJesAuthentication with LazyLogging {

  override def stdoutStderr: CallLogs = backend.stdoutStderr(jobDescriptor)
}
