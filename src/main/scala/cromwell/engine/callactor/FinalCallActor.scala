package cromwell.engine.callactor

import akka.actor.ActorRef
import cromwell.engine.WorkflowDescriptor
import cromwell.engine.backend.{FinalCallJobDescriptor, BackendCallJobDescriptor}
import cromwell.engine.callexecution.CallExecutionActor
import cromwell.engine.workflow.FinalCallKey

case class FinalCallActor(jobDescriptor: FinalCallJobDescriptor) extends CallActor[FinalCallJobDescriptor] {
  override protected lazy val callExecutionActor: ActorRef = context.actorOf(CallExecutionActor.props(jobDescriptor.key.scope))
}
