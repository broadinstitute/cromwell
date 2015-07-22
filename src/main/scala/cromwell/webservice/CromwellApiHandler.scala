package cromwell.webservice

import akka.actor.{Actor, ActorRef, Props}
import akka.event.Logging
import akka.pattern.ask
import akka.util.Timeout
import cromwell.binding.{WdlJson, WdlSource, WorkflowRawInputs}
import cromwell.engine._
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowManagerActorMessage
import cromwell.parser.WdlParser.SyntaxError
import cromwell.webservice.CromwellApiHandler._
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.{binding, engine}
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}

object CromwellApiHandler {

  def props(workflowManagerActorRef: ActorRef): Props = {
    Props(new CromwellApiHandler(workflowManagerActorRef))
  }

  sealed trait WorkflowManagerMessage

  case class SubmitWorkflow(wdlSource: WdlSource, wdlJson: WdlJson, inputs: WorkflowRawInputs) extends WorkflowManagerMessage

  case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerMessage

  case class WorkflowOutputs(id: WorkflowId) extends WorkflowManagerMessage

  case class WorkflowAbort(id: WorkflowId) extends WorkflowManagerMessage
  
  // Builds an internal message for the workflow manager, which will then respond to an ask with an
  // Option[WorkflowState].
  type InternalMessageBuilder = WorkflowId => WorkflowManagerActorMessage
}

class CromwellApiHandler(workflowManager: ActorRef) extends Actor {

  import WorkflowJsonSupport._
  import context.dispatcher

  implicit val timeout = Timeout(2 seconds)

  val log = Logging(context.system, classOf[CromwellApiHandler])

  override def receive = {

    case WorkflowStatus(id) =>
      val futureStatus = queryWorkflowManager(WorkflowManagerActor.WorkflowStatus(id))
      futureStatus onSuccess {
        case None => context.parent ! RequestComplete(StatusCodes.NotFound)
        case Some(workflowState) =>
          context.parent ! RequestComplete(StatusCodes.OK, WorkflowStatusResponse(id.toString, workflowState.toString))
      }

    case WorkflowAbort(id) =>
      val futureStatus = queryWorkflowManager(WorkflowManagerActor.WorkflowStatus(id))
      futureStatus onSuccess {
        case None =>
          context.parent ! RequestComplete(StatusCodes.NotFound, s"Workflow ID '$id' not found.")
        case Some(workflowState) if workflowState.isTerminal =>
          val message = s"Workflow ID '$id' is in terminal state '$workflowState' and cannot be aborted."
          context.parent ! RequestComplete(StatusCodes.Forbidden, message)
        case Some(workflowState) =>
          val futureAbortedStatus = queryWorkflowManager(WorkflowManagerActor.WorkflowAbort(id))
          futureAbortedStatus onComplete {
            case Failure(e) =>
              log.error(e, e.getMessage)
              context.parent ! RequestComplete(StatusCodes.InternalServerError, e.getMessage)
            case Success(status) =>
              context.parent ! RequestComplete(StatusCodes.OK, WorkflowAbortResponse(id.toString, WorkflowAborted.toString))
          }
      }


    case SubmitWorkflow(wdlSource, wdlJson, inputs) =>
      val workflowManagerResponseFuture = ask(workflowManager, WorkflowManagerActor.SubmitWorkflow(wdlSource, wdlJson, inputs)).mapTo[WorkflowId]
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
      val eventualWorkflowOutputs = ask(workflowManager, WorkflowManagerActor.WorkflowOutputs(id)).mapTo[binding.WorkflowOutputs]
      eventualWorkflowOutputs onComplete {
        case Success(outputs) => context.parent ! RequestComplete(StatusCodes.OK, WorkflowOutputResponse(id.toString, outputs))
        case Failure(ex: WorkflowManagerActor.WorkflowNotFoundException) => context.parent ! RequestComplete(StatusCodes.NotFound)
        case Failure(ex) => context.parent ! RequestComplete(StatusCodes.InternalServerError, ex.getMessage)
      }
  }

  /**
   * Query the workflow manager with the specified message and return a `Future[Option[WorkflowState]]`.  Add
   * an `onFailure` handler for the `Future` to message the `context.parent` with a `StatusCodes.InternalServerError`
   * if the `Future` fails.
   *
   */
  private def queryWorkflowManager(message: WorkflowManagerActorMessage): Future[Option[WorkflowState]] = {
    val future = ask(workflowManager, message).mapTo[Option[WorkflowState]]
    future onFailure {
      case ex =>
        log.error(ex, ex.getMessage)
        context.parent ! RequestComplete(StatusCodes.InternalServerError, ex.getMessage)
    }
    future
  }
}
