package cromwell.engine.callactor

import cromwell.engine.{AbortFunction, AbortRegistrationFunction}
import cromwell.engine.backend.OldStyleBackendCallJobDescriptor
import cromwell.engine.callexecution.OldStyleCallExecutionActor

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class OldStyleBackendCallActor(jobDescriptor: OldStyleBackendCallJobDescriptor) extends OldStyleCallActor[OldStyleBackendCallJobDescriptor] {

  override lazy val callExecutionActor = {
    def registerAbortFunction(abortFunction: AbortFunction): Unit = {
      self ! OldStyleCallActor.RegisterCallAbortFunction(abortFunction)
    }
    val backendCall = jobDescriptor.copy(abortRegistrationFunction = Option(AbortRegistrationFunction(registerAbortFunction)))
    val executionActorName = s"CallExecutionActor-${workflowDescriptor.id}-${call.unqualifiedName}"
    context.actorOf(OldStyleCallExecutionActor.props(backendCall), executionActorName)
  }
}
