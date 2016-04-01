package cromwell.engine.callactor

import akka.actor.ActorRef
import cromwell.engine.backend.FinalCallJobDescriptor
import cromwell.engine.callexecution.CallExecutionActor

case class FinalCallActor(jobDescriptor: FinalCallJobDescriptor) extends CallActor[FinalCallJobDescriptor] {
  override protected lazy val callExecutionActor: ActorRef = context.actorOf(
    CallExecutionActor.props(jobDescriptor.call, jobDescriptor.workflowMetadataResponse))
}
