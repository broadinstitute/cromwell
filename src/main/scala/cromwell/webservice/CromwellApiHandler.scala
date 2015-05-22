package cromwell.webservice

import akka.actor.{Actor, ActorRef, Props}
import java.util.UUID

import akka.actor.{ActorRef, Actor, Props}
import akka.pattern.ask
import akka.util.Timeout
import cromwell.binding
import cromwell.binding._
import cromwell.engine._
import cromwell.webservice.CromwellApiHandler.{WorkflowOutputs, WorkflowStatus}
import cromwell.webservice.PerRequest.RequestComplete
import spray.http.StatusCodes

import scala.concurrent.duration._
import scala.util.{Failure, Success}

object CromwellApiHandler {
  sealed trait WorkflowManagerMessage
  case class SubmitWorkflow(wdl: WdlSource, inputs: WorkflowRawInputs) extends WorkflowManagerMessage
  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerMessage
  case class WorkflowOutputs(id: WorkflowId) extends WorkflowManagerMessage

  def props(workflowManagerActorRef : ActorRef): Props = {
    Props(new CromwellApiHandler(workflowManagerActorRef))
  }
}

class CromwellApiHandler(workflowManager : ActorRef) extends Actor {
  import context.dispatcher
  implicit val timeout = Timeout(2.seconds)

  override def receive = {
    case WorkflowStatus(id) => handleWorkflowStatus(id)
    case WorkflowOutputs(id) => handleWorkflowOutputs(id)
  }

  def handleWorkflowStatus(id: UUID): Unit = {
    val eventualWorkflowStatus = ask(workflowManager, WorkflowManagerActor.WorkflowStatus(id)).mapTo[Option[WorkflowState]]

    eventualWorkflowStatus onComplete {
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

  def handleWorkflowOutputs(id: UUID): Unit = {
    // FIXME: Abstract out similar code from handleWorkflowStatus - perhaps a func to take the Future and a func to process for Success?
    // FIXME: If the general pattern is Future[Option[T]] I think it should be a func that takes that and a fun to proess Success(Some)) case
    val eventualWorkflowOutputs = ask(workflowManager, WorkflowManagerActor.WorkflowOutputs(id)).mapTo[Option[binding.WorkflowOutputs]]
    eventualWorkflowOutputs onComplete {
      case Success(outputs) => outputs match {
        case Some(x) =>
          val outputMap = x mapValues {_.toString}
          context.parent ! RequestComplete(StatusCodes.OK, WorkflowOutputResponse(id.toString, outputMap))
        case None => context.parent ! RequestComplete(StatusCodes.NotFound, None)
      }
      case Failure(ex) => context.parent ! RequestComplete(StatusCodes.InternalServerError, ex.getMessage)
    }
  }
}
