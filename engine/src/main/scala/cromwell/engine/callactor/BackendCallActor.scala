package cromwell.engine.callactor

import cromwell.engine.backend.BackendCallJobDescriptor
import cromwell.engine.callexecution.CallExecutionActor
import cromwell.engine.{AbortFunction, AbortRegistrationFunction}

case class BackendCallActor(jobDescriptor: BackendCallJobDescriptor) extends CallActor {

  override lazy val callExecutionActor = {
    def registerAbortFunction(abortFunction: AbortFunction): Unit = {
      self ! CallActor.RegisterCallAbortFunction(abortFunction)
    }
    val backendCall = jobDescriptor.copy(abortRegistrationFunction = Option(AbortRegistrationFunction(registerAbortFunction)))
    val executionActorName = s"CallExecutionActor-${workflowDescriptor.id}-${call.unqualifiedName}"
    context.actorOf(CallExecutionActor.props(backendCall), executionActorName)
  }
}
