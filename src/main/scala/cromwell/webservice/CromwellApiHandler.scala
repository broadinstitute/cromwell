package cromwell.webservice

import akka.actor.{Actor, ActorRef, Props}
import akka.pattern.ask
import akka.util.Timeout
import cromwell.binding.{WdlSource, WorkflowRawInputs}
import cromwell.engine._
import cromwell.parser.WdlParser.SyntaxError
import cromwell.webservice.CromwellApiHandler._
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.{binding, engine}
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._

import scala.concurrent.duration._
import scala.util.{Failure, Success}

object CromwellApiHandler {

  def props(workflowManagerActorRef: ActorRef): Props = {
    Props(new CromwellApiHandler(workflowManagerActorRef))
  }

  sealed trait WorkflowManagerMessage

  case class SubmitWorkflow(wdl: WdlSource, inputs: WorkflowRawInputs) extends WorkflowManagerMessage

  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerMessage

  case class WorkflowOutputs(id: WorkflowId) extends WorkflowManagerMessage
}

class CromwellApiHandler(workflowManager: ActorRef) extends Actor {

  import WorkflowJsonSupport._
  import cromwell.binding.values.WdlValueJsonFormatter._
  import context.dispatcher

  implicit val timeout = Timeout(2.seconds)

  override def receive = {
    case WorkflowStatus(id) =>
      val workflowManagerResponseFuture = ask(workflowManager, WorkflowManagerActor.WorkflowStatus(id)).mapTo[Option[WorkflowState]]

      workflowManagerResponseFuture.onComplete {
        case Success(state) =>
          state match {
            case Some(x) =>
              context.parent ! RequestComplete(StatusCodes.OK, WorkflowStatusResponse(id.toString, x.toString))
            case None => context.parent ! RequestComplete(StatusCodes.NotFound)
          }

        case Failure(ex) =>
          context.parent ! RequestComplete(StatusCodes.InternalServerError, ex.getMessage)
      }

    case SubmitWorkflow(wdl, inputs) =>
      val workflowManagerResponseFuture = ask(workflowManager, WorkflowManagerActor.SubmitWorkflow(wdl, inputs)).mapTo[WorkflowId]
      workflowManagerResponseFuture.onComplete {
        case Success(id) =>
          context.parent ! RequestComplete(StatusCodes.Created, WorkflowSubmitResponse(id.toString, engine.WorkflowSubmitted.toString))

        case Failure(ex) =>
          ex match {
            case _: SyntaxError =>
              context.parent ! RequestComplete(StatusCodes.BadRequest, ex.getMessage)
            case _ =>
              context.parent ! RequestComplete(StatusCodes.InternalServerError, ex.getMessage)
          }
      }

    case WorkflowOutputs(id) =>
      val eventualWorkflowOutputs = ask(workflowManager, WorkflowManagerActor.WorkflowOutputs(id)).mapTo[Option[binding.WorkflowOutputs]]
      eventualWorkflowOutputs onComplete {
        case Success(outputs) => outputs match {
          case Some(x) =>
            context.parent ! RequestComplete(StatusCodes.OK, WorkflowOutputResponse(id.toString, x))
          case None => context.parent ! RequestComplete(StatusCodes.NotFound)
        }
        case Failure(ex) => context.parent ! RequestComplete(StatusCodes.InternalServerError, ex.getMessage)
      }
  }
}
