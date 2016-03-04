package cromwell.engine.callactor

import cromwell.engine.AbortRegistrationFunction
import cromwell.engine.backend.BackendCallJobDescriptor
import cromwell.engine.callexecution.CallExecutionActor
import cromwell.engine.workflow.BackendCallKey

case class BackendCallActor(jobDescriptor: BackendCallJobDescriptor) extends CallActor[BackendCallJobDescriptor] {

  override lazy val callExecutionActor = {
    val backendCall = jobDescriptor.backend.bindCall(jobDescriptor, Option(AbortRegistrationFunction(registerAbortFunction)))
    val executionActorName = s"CallExecutionActor-${workflowDescriptor.id}-${call.unqualifiedName}"
    context.actorOf(CallExecutionActor.props(backendCall), executionActorName)
  }
}
