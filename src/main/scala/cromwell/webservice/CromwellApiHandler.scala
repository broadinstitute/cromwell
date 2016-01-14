package cromwell.webservice

import akka.actor.{Actor, ActorRef, Props}
import akka.event.Logging
import akka.pattern.ask
import akka.util.Timeout
import cromwell.engine._
import cromwell.engine.backend.{Backend, CallLogs}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowManagerActorMessage
import cromwell.webservice.CromwellApiHandler.{CallOutputs, WorkflowOutputs, _}
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.engine
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._
import spray.json._

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}

object CromwellApiHandler {
  def props(workflowManagerActorRef: ActorRef): Props = {
    Props(new CromwellApiHandler(workflowManagerActorRef))
  }

  sealed trait WorkflowManagerMessage

  final case class WorkflowSubmit(source: WorkflowSourceFiles) extends WorkflowManagerMessage
  final case class WorkflowQuery(parameters: Seq[(String, String)]) extends WorkflowManagerMessage
  final case class WorkflowStatus(id: WorkflowId) extends WorkflowManagerMessage
  final case class WorkflowOutputs(id: WorkflowId) extends WorkflowManagerMessage
  final case class WorkflowAbort(id: WorkflowId) extends WorkflowManagerMessage
  final case class CallOutputs(id: WorkflowId, callFqn: String) extends WorkflowManagerMessage
  final case class CallStdoutStderr(id: WorkflowId, callFqn: String) extends WorkflowManagerMessage
  final case class WorkflowStdoutStderr(id: WorkflowId) extends WorkflowManagerMessage
  final case class CallCaching(id: WorkflowId, parameters: QueryParameters, callName: Option[String]) extends WorkflowManagerMessage
  final case class WorkflowMetadata(id: WorkflowId) extends WorkflowManagerMessage
}

class CromwellApiHandler(workflowManager: ActorRef) extends Actor {
  import WorkflowJsonSupport._
  import context.dispatcher

  implicit val timeout = Timeout(2 seconds)
  val log = Logging(context.system, classOf[CromwellApiHandler])

  def workflowNotFound(id: WorkflowId) = RequestComplete(StatusCodes.NotFound, APIResponse.error(new Throwable(s"Workflow '$id' not found.")))
  def callNotFound(callFqn: String, id: WorkflowId) = {
    RequestComplete(StatusCodes.NotFound, APIResponse.error(new Throwable(s"Call $callFqn not found for workflow '$id'.")))
  }

  override def receive = {
    case WorkflowStatus(id) =>
      val futureStatus = queryWorkflowManager(WorkflowManagerActor.WorkflowStatus(id))
      futureStatus onSuccess {
        case None => context.parent ! RequestComplete(StatusCodes.NotFound)
        case Some(workflowState) =>
          context.parent ! RequestComplete(StatusCodes.OK, WorkflowStatusResponse(id.toString, workflowState.toString))
      }
    case WorkflowQuery(parameters) =>
      val futureQuery = ask(workflowManager, WorkflowManagerActor.WorkflowQuery(parameters)).mapTo[WorkflowQueryResponse]
      futureQuery onComplete {
        case Success(response) => context.parent ! RequestComplete(StatusCodes.OK, response)
        case Failure(ex) =>
          ex match {
            case _: IllegalArgumentException =>
              context.parent ! RequestComplete(StatusCodes.BadRequest, APIResponse.fail(ex))
            case _ =>
              context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(ex))
          }
      }
    case WorkflowAbort(id) =>
      val futureStatus = queryWorkflowManager(WorkflowManagerActor.WorkflowStatus(id))
      futureStatus onSuccess {
        case None =>
          context.parent ! workflowNotFound(id)
        case Some(workflowState) if workflowState.isTerminal =>
          val message = s"Workflow ID '$id' is in terminal state '$workflowState' and cannot be aborted."
          context.parent ! RequestComplete(StatusCodes.Forbidden, APIResponse.fail(new Throwable(message)))
        case Some(workflowState) =>
          val futureAbortedStatus = queryWorkflowManager(WorkflowManagerActor.WorkflowAbort(id))
          futureAbortedStatus onComplete {
            case Failure(e) =>
              log.error(e, e.getMessage)
              context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(e))
            case Success(status) =>
              context.parent ! RequestComplete(StatusCodes.OK, WorkflowAbortResponse(id.toString, WorkflowAborted.toString))
          }
      }
    case WorkflowSubmit(source) =>
      val workflowManagerResponseFuture = ask(workflowManager, WorkflowManagerActor.SubmitWorkflow(source)).mapTo[WorkflowId]
      workflowManagerResponseFuture.onComplete {
        case Success(id) =>
          context.parent ! RequestComplete(StatusCodes.Created, WorkflowSubmitResponse(id.toString, engine.WorkflowSubmitted.toString))
        case Failure(ex) =>
          ex match {
            case _: IllegalArgumentException => context.parent ! RequestComplete(StatusCodes.BadRequest, APIResponse.fail(ex))
            case _ => context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(ex))
          }
      }
    case WorkflowOutputs(id) =>
      val eventualWorkflowOutputs = ask(workflowManager, WorkflowManagerActor.WorkflowOutputs(id)).mapTo[engine.WorkflowOutputs]
      eventualWorkflowOutputs onComplete {
        case Success(outputs) => context.parent ! RequestComplete(StatusCodes.OK, WorkflowOutputResponse(id.toString, outputs.mapToValues))
        case Failure(ex: WorkflowManagerActor.WorkflowNotFoundException) => context.parent ! workflowNotFound(id)
        case Failure(ex) => context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(ex))
      }
    case CallOutputs(id, callFqn) =>
      val eventualCallOutputs = ask(workflowManager, WorkflowManagerActor.CallOutputs(id, callFqn)).mapTo[engine.CallOutputs]
      eventualCallOutputs onComplete {
        case Success(outputs) if outputs.nonEmpty => context.parent ! RequestComplete(StatusCodes.OK, CallOutputResponse(id.toString, callFqn, outputs.mapToValues))
        case Failure(ex: WorkflowManagerActor.WorkflowNotFoundException) => context.parent ! workflowNotFound(id)
        case Failure(ex: WorkflowManagerActor.CallNotFoundException) => context.parent ! callNotFound(callFqn, id)
        case Failure(ex) => context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(ex))
      }
    case CallStdoutStderr(id, callFqn) =>
      val eventualCallLogs = ask(workflowManager, WorkflowManagerActor.CallStdoutStderr(id, callFqn)).mapTo[Seq[CallLogs]]
      eventualCallLogs onComplete {
        case Success(logs) => context.parent ! RequestComplete(StatusCodes.OK, CallStdoutStderrResponse(id.toString, Map(callFqn -> logs)))
        case Failure(ex: WorkflowManagerActor.WorkflowNotFoundException) => context.parent ! workflowNotFound(id)
        case Failure(ex: WorkflowManagerActor.CallNotFoundException) => context.parent ! callNotFound(callFqn, id)
        case Failure(ex: Backend.StdoutStderrException) => context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(ex))
        case Failure(ex) => context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(ex))
      }
    case WorkflowStdoutStderr(id) =>
      val eventualCallLogs = ask(workflowManager, WorkflowManagerActor.WorkflowStdoutStderr(id)).mapTo[Map[String, Seq[CallLogs]]]
      eventualCallLogs onComplete {
        case Success(logs) =>
          context.parent ! RequestComplete(StatusCodes.OK, CallStdoutStderrResponse(id.toString, logs))
        case Failure(ex: WorkflowManagerActor.WorkflowNotFoundException) => context.parent ! workflowNotFound(id)
        case Failure(ex) => context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(ex))
      }
    case WorkflowMetadata(id) =>
      val eventualMetadataResponse = ask(workflowManager, WorkflowManagerActor.WorkflowMetadata(id)).mapTo[WorkflowMetadataResponse]
      eventualMetadataResponse onComplete {
        case Success(metadata) => context.parent ! RequestComplete(StatusCodes.OK, metadata)
        case Failure(ex: WorkflowManagerActor.WorkflowNotFoundException) => context.parent ! workflowNotFound(id)
        case Failure(ex) => context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(ex))
      }
    case CallCaching(id, parameters, callName) =>
      val eventualUpdateCount = ask(workflowManager, WorkflowManagerActor.CallCaching(id, parameters, callName)).mapTo[Int]
      eventualUpdateCount onComplete {
        case Success(updateCount) => context.parent ! RequestComplete(StatusCodes.OK, CallCachingResponse(updateCount))
        case Failure(ex: WorkflowManagerActor.WorkflowNotFoundException) => context.parent ! workflowNotFound(id)
        case Failure(ex: IllegalArgumentException) => context.parent ! RequestComplete(StatusCodes.BadRequest, APIResponse.fail(ex))
        case Failure(ex) => context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(ex))
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
        context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(ex))
    }

    future
  }
}
