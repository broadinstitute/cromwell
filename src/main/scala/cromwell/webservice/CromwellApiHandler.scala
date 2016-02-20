package cromwell.webservice

import akka.actor.{Actor, ActorRef, Props}
import akka.event.Logging
import akka.util.Timeout
import cromwell.engine
import cromwell.engine._
import cromwell.engine.backend.CallLogs
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.{CallNotFoundException, WorkflowNotFoundException}
import cromwell.webservice.PerRequest.RequestComplete
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._

import scala.concurrent.duration._
import scala.language.postfixOps

object CromwellApiHandler {
  def props(workflowManagerActorRef: ActorRef): Props = {
    Props(new CromwellApiHandler(workflowManagerActorRef))
  }

  sealed trait ApiHandlerMessage

  final case class ApiHandlerWorkflowSubmit(source: WorkflowSourceFiles) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowQuery(parameters: Seq[(String, String)]) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowStatus(id: WorkflowId) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowOutputs(id: WorkflowId) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowAbort(id: WorkflowId) extends ApiHandlerMessage
  final case class ApiHandlerCallOutputs(id: WorkflowId, callFqn: String) extends ApiHandlerMessage
  final case class ApiHandlerCallStdoutStderr(id: WorkflowId, callFqn: String) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowStdoutStderr(id: WorkflowId) extends ApiHandlerMessage
  final case class ApiHandlerCallCaching(id: WorkflowId, parameters: QueryParameters, callName: Option[String]) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowMetadata(id: WorkflowId) extends ApiHandlerMessage

  sealed trait WorkflowManagerResponse

  sealed trait WorkflowManagerSuccessResponse extends WorkflowManagerResponse

  sealed trait WorkflowManagerFailureResponse extends WorkflowManagerResponse {
    def failure: Throwable
  }

  final case class WorkflowManagerSubmitSuccess(id: WorkflowId) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerSubmitFailure(override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerWorkflowOutputsSuccess(id: WorkflowId, outputs: engine.WorkflowOutputs) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerWorkflowOutputsFailure(id: WorkflowId, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerStatusSuccess(id: WorkflowId, state: WorkflowState) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerStatusFailure(id: WorkflowId, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerAbortSuccess(id: WorkflowId) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerAbortFailure(id: WorkflowId, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerQuerySuccess(response: WorkflowQueryResponse) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerQueryFailure(override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerCallOutputsSuccess(id: WorkflowId, callFqn: FullyQualifiedName, outputs: engine.CallOutputs) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerCallOutputsFailure(id: WorkflowId, callFqn: FullyQualifiedName, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerCallStdoutStderrSuccess(id: WorkflowId, callFqn: FullyQualifiedName, logs: Seq[CallLogs]) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerCallStdoutStderrFailure(id: WorkflowId, callFqn: FullyQualifiedName, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerWorkflowStdoutStderrSuccess(id: WorkflowId, logs: Map[FullyQualifiedName, Seq[CallLogs]]) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerWorkflowStdoutStderrFailure(id: WorkflowId, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerWorkflowMetadataSuccess(id: WorkflowId, response: WorkflowMetadataResponse) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerWorkflowMetadataFailure(id: WorkflowId, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerCallCachingSuccess(id: WorkflowId, updateCount: Int) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerCallCachingFailure(id: WorkflowId, override val failure: Throwable) extends WorkflowManagerFailureResponse
}

class CromwellApiHandler(workflowManager: ActorRef) extends Actor {
  import CromwellApiHandler._
  import WorkflowJsonSupport._

  implicit val timeout = Timeout(2 seconds)
  val log = Logging(context.system, classOf[CromwellApiHandler])

  def workflowNotFound(id: WorkflowId) = RequestComplete(StatusCodes.NotFound, APIResponse.error(new Throwable(s"Workflow '$id' not found.")))
  def callNotFound(callFqn: String, id: WorkflowId) = {
    RequestComplete(StatusCodes.NotFound, APIResponse.error(new Throwable(s"Call $callFqn not found for workflow '$id'.")))
  }

  private def error(t: Throwable)(f: Throwable => RequestComplete[_]): Unit = context.parent ! f(t)

  override def receive = {
    case ApiHandlerWorkflowStatus(id) => workflowManager ! WorkflowManagerActor.WorkflowStatus(id)
    case WorkflowManagerStatusSuccess(id, state) => context.parent ! RequestComplete(StatusCodes.OK, WorkflowStatusResponse(id.toString, state.toString))
    case WorkflowManagerStatusFailure(_, e) =>
      error(e) {
        case _: WorkflowNotFoundException => RequestComplete(StatusCodes.NotFound, e)
        case _ => RequestComplete(StatusCodes.InternalServerError, e)
      }

    case ApiHandlerWorkflowQuery(parameters) => workflowManager ! WorkflowManagerActor.WorkflowQuery(parameters)
    case WorkflowManagerQuerySuccess(response) => context.parent ! RequestComplete(StatusCodes.OK, response)
    case WorkflowManagerQueryFailure(e) =>
      error(e) {
        case _: IllegalArgumentException => RequestComplete(StatusCodes.BadRequest, APIResponse.fail(e))
        case _ => RequestComplete(StatusCodes.InternalServerError, APIResponse.error(e))
      }

    case ApiHandlerWorkflowAbort(id) => workflowManager ! WorkflowManagerActor.WorkflowAbort(id)
    case WorkflowManagerAbortSuccess(id) =>
      context.parent ! RequestComplete(StatusCodes.OK, WorkflowAbortResponse(id.toString, WorkflowAborted.toString))
    case WorkflowManagerAbortFailure(_, e) =>
      error(e) {
        case _: WorkflowNotFoundException => RequestComplete(StatusCodes.NotFound, APIResponse.error(e))
        case _: IllegalStateException => RequestComplete(StatusCodes.Forbidden, APIResponse.error(e))
        case _ => RequestComplete(StatusCodes.InternalServerError, APIResponse.error(e))
      }

    case ApiHandlerWorkflowSubmit(source) => workflowManager ! WorkflowManagerActor.ValidateAndSubmitWorkflow(source)
    case WorkflowManagerSubmitSuccess(id) =>
      context.parent ! RequestComplete(StatusCodes.Created, WorkflowSubmitResponse(id.toString, engine.WorkflowSubmitted.toString))
    case WorkflowManagerSubmitFailure(exception) =>
      exception match {
        case reason: IllegalStateException => context.parent ! RequestComplete(StatusCodes.BadRequest, APIResponse.fail(exception))
        case _ => context.parent ! RequestComplete(StatusCodes.InternalServerError, APIResponse.error(exception))
      }
    case ApiHandlerWorkflowOutputs(id) => workflowManager ! WorkflowManagerActor.WorkflowOutputs(id)
    case WorkflowManagerWorkflowOutputsSuccess(id, outputs) =>
      context.parent ! RequestComplete(StatusCodes.OK, WorkflowOutputResponse(id.toString, outputs.mapToValues))
    case WorkflowManagerWorkflowOutputsFailure(id, e) =>
      error(e) {
        case _: WorkflowNotFoundException => workflowNotFound(id)
        case _ => RequestComplete(StatusCodes.InternalServerError, APIResponse.error(e))
      }

    case ApiHandlerCallOutputs(id, callFqn) => workflowManager ! WorkflowManagerActor.CallOutputs(id, callFqn)
    case WorkflowManagerCallOutputsSuccess(id, callFqn, outputs) =>
      context.parent ! RequestComplete(StatusCodes.OK, CallOutputResponse(id.toString, callFqn, outputs.mapToValues))
    case WorkflowManagerCallOutputsFailure(id, callFqn, e) =>
      error(e) {
        case _: WorkflowNotFoundException => workflowNotFound(id)
        case _: CallNotFoundException => callNotFound(callFqn, id)
        case _ => RequestComplete(StatusCodes.InternalServerError, APIResponse.error(e))
      }

    case ApiHandlerCallStdoutStderr(id, callFqn) => workflowManager ! WorkflowManagerActor.CallStdoutStderr(id, callFqn)
    case WorkflowManagerCallStdoutStderrSuccess(id, callFqn, logs) =>
      context.parent ! RequestComplete(StatusCodes.OK, CallStdoutStderrResponse(id.toString, Map(callFqn -> logs)))
    case WorkflowManagerCallStdoutStderrFailure(id, callFqn, e) =>
      error(e) {
        case _: WorkflowNotFoundException => workflowNotFound(id)
        case _: CallNotFoundException => callNotFound(callFqn, id)
          //TODO: Pending on Issue #430
//        case _: StdoutStderrException =>
//          RequestComplete(StatusCodes.InternalServerError, APIResponse.error(e)) //"Backend.StdoutStderrException => RequestComplete(StatusCodes.InternalServerError, APIResponse.error(e))"
        case _ => RequestComplete(StatusCodes.InternalServerError, APIResponse.error(e))
      }

    case ApiHandlerWorkflowStdoutStderr(id) => workflowManager ! WorkflowManagerActor.WorkflowStdoutStderr(id)
    case WorkflowManagerWorkflowStdoutStderrSuccess(id, callLogs) =>
          context.parent ! RequestComplete(StatusCodes.OK, CallStdoutStderrResponse(id.toString, callLogs))
    case WorkflowManagerWorkflowStdoutStderrFailure(id, e) =>
      error(e) {
        case _: WorkflowNotFoundException => workflowNotFound(id)
        case _ => RequestComplete(StatusCodes.InternalServerError, APIResponse.error(e))
      }

    case ApiHandlerWorkflowMetadata(id) => workflowManager ! WorkflowManagerActor.WorkflowMetadata(id)
    case WorkflowManagerWorkflowMetadataSuccess(id, response) => context.parent ! RequestComplete(StatusCodes.OK, response)
    case WorkflowManagerWorkflowMetadataFailure(id, e) =>
      error(e) {
        case _: WorkflowNotFoundException => workflowNotFound(id)
        case _ => RequestComplete(StatusCodes.InternalServerError, APIResponse.error(e))
      }

    case ApiHandlerCallCaching(id, parameters, callName) => workflowManager ! WorkflowManagerActor.CallCaching(id, parameters, callName)
    case WorkflowManagerCallCachingSuccess(id, updateCount) => context.parent ! RequestComplete(StatusCodes.OK, CallCachingResponse(updateCount))
    case WorkflowManagerCallCachingFailure(id, e) =>
      error(e) {
        case _: WorkflowNotFoundException => workflowNotFound(id)
        case _: IllegalArgumentException => RequestComplete(StatusCodes.BadRequest, APIResponse.fail(e))
        case _ => RequestComplete(StatusCodes.InternalServerError, APIResponse.error(e))
      }
  }
}
