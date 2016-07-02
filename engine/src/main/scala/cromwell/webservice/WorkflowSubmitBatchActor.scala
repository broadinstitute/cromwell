package cromwell.webservice

import akka.actor.{Actor, ActorRef, Props}
import cromwell.core.WorkflowSourceFiles
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.webservice.CromwellApiHandler.{WorkflowManagerBatchSubmitResponse, WorkflowManagerResponse}

/** Submits a sequence of sources, then messages back the respective sequence of responses. */
class WorkflowSubmitBatchActor(responseActor: ActorRef, requestHandlerActor: ActorRef,
                               sources: Seq[WorkflowSourceFiles]) extends Actor {

  private val responses = Array.fill[Option[WorkflowManagerResponse]](sources.size)(None)

  override def preStart() = {
    sources.zipWithIndex foreach {
      case (source, index) =>
        context.actorOf(
          Props(new WorkflowSubmitIndexedActor(self, requestHandlerActor, source, index)),
          "WorkflowSubmitIndexedActor_" + index)
    }
  }

  override def receive = {
    case WorkflowManagerBatchIndexedResponse(response, index) =>
      responses(index) = Option(response)
      if (responses.forall(_.isDefined)) {
        responseActor ! WorkflowManagerBatchSubmitResponse(responses.map(_.get))
        context stop self
      }
  }
}

case class WorkflowManagerBatchIndexedResponse(response: WorkflowManagerResponse, index: Int)

class WorkflowSubmitIndexedActor(batchSubmitActor: ActorRef, requestHandlerActor: ActorRef,
                                 source: WorkflowSourceFiles, index: Int) extends Actor {

  override def preStart() = {
    requestHandlerActor ! WorkflowManagerActor.SubmitWorkflowCommand(source)
  }

  override def receive = {
    case response: WorkflowManagerResponse =>
      batchSubmitActor ! WorkflowManagerBatchIndexedResponse(response, index)
      context stop self
  }
}
