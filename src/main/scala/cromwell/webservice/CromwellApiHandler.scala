package cromwell.webservice

import akka.actor.{Actor, ActorRef, Props}
import akka.pattern.ask
import akka.util.Timeout
import cromwell.binding.{WdlSource, WorkflowRawInputs}
import cromwell.engine
import cromwell.engine._
import cromwell.parser.WdlParser.SyntaxError
import cromwell.webservice.CromwellApiHandler._
import cromwell.webservice.PerRequest.RequestComplete
import spray.http.StatusCodes

import scala.concurrent.duration._
import scala.util.{Failure, Success}

object CromwellApiHandler {
  sealed trait WorkflowManagerMessage
  case class SubmitWorkflow(wdl: WdlSource, inputs: WorkflowRawInputs) extends WorkflowManagerMessage
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

    case SubmitWorkflow(wdl, inputs) =>
      implicit val timeout = Timeout(2.seconds)
      val workflowManagerResponseFuture = ask(workflowManagerActorRef, WorkflowManagerActor.SubmitWorkflow(wdl, inputs)).mapTo[WorkflowId]

      workflowManagerResponseFuture.onComplete {
        case Success(id) =>
          context.parent ! RequestComplete (StatusCodes.Created, WorkflowSubmitResponse(id.toString, engine.WorkflowSubmitted.toString))

        case Failure(ex) =>
          ex match {
            case _:SyntaxError =>
              context.parent ! RequestComplete(StatusCodes.BadRequest, ex.getMessage)
            case _ =>
              context.parent ! RequestComplete(StatusCodes.InternalServerError, ex.getMessage)
          }
      }

  }
}
