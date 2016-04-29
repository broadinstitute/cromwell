package cromwell.engine.callactor

import akka.actor.ActorRef
import cromwell.engine.backend.FinalCallJobDescriptor
import cromwell.engine.callexecution.OldStyleCallExecutionActor

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class OldStyleFinalCallActor(jobDescriptor: FinalCallJobDescriptor) extends OldStyleCallActor[FinalCallJobDescriptor] {
  override protected lazy val callExecutionActor: ActorRef = context.actorOf(
    OldStyleCallExecutionActor.props(jobDescriptor.call, jobDescriptor.workflowMetadataResponse))
}
