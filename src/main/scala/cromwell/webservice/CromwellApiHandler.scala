package cromwell.webservice

import akka.actor.{ActorRef, Actor, Props}
import akka.pattern.ask
import akka.util.Timeout
import cromwell.binding._
import cromwell.engine._
import cromwell.webservice.CromwellApiHandler.WorkflowStatus
import cromwell.webservice.PerRequest.RequestComplete
import spray.http.StatusCodes

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.util.{Failure, Success}

object CromwellApiHandler {
  sealed trait WorkflowManagerMessage
  case class SubmitWorkflow(wdl: WdlSource, inputs: WorkflowInputs) extends WorkflowManagerMessage
  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerMessage

  def props(workflowManagerActorRef : ActorRef): Props = {
    Props(new CromwellApiHandler(workflowManagerActorRef))
  }
}

class CromwellApiHandler(workflowManagerActorRef : ActorRef) extends Actor {
  import context.dispatcher

  override def receive = {
    case WorkflowStatus(id) =>
      implicit val timeout = Timeout(2.seconds)
      val workflowManagerResponseFuture = ask(workflowManagerActorRef, WorkflowManagerActor.WorkflowStatus(id)).mapTo[Option[WorkflowState]]

      workflowManagerResponseFuture.onComplete {
        case Success(state) =>
          state match {
            case Some(x) =>
              context.parent ! RequestComplete(StatusCodes.OK, WorkflowStatusResponse(id.toString, x.toString))
            case None => context.parent ! RequestComplete(StatusCodes.NotFound, None)
          }

        case Failure(ex) =>
          context.parent ! RequestComplete(StatusCodes.InternalServerError, ex.getMessage)
      }
  }
}
