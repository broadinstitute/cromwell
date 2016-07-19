package cromwell.webservice

import akka.actor.{Actor, ActorRef, Props}
import akka.event.Logging
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import cromwell.core._
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.services.MetadataServiceActor.{QueryMetadata, WorkflowQueryResponse}
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.core
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor
import spray.http.{StatusCodes, Uri}
import spray.httpx.SprayJsonSupport._

import scala.concurrent.duration._
import scala.language.postfixOps
import scalaz.NonEmptyList

object CromwellApiHandler {
  def props(requestHandlerActor: ActorRef): Props = {
    Props(new CromwellApiHandler(requestHandlerActor))
  }

  sealed trait ApiHandlerMessage

  final case class ApiHandlerWorkflowSubmit(source: WorkflowSourceFiles) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowSubmitBatch(sources: NonEmptyList[WorkflowSourceFiles]) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowQuery(uri: Uri, parameters: Seq[(String, String)]) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowStatus(id: WorkflowId) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowOutputs(id: WorkflowId) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowAbort(id: WorkflowId) extends ApiHandlerMessage
  final case class ApiHandlerCallOutputs(id: WorkflowId, callFqn: String) extends ApiHandlerMessage
  final case class ApiHandlerCallStdoutStderr(id: WorkflowId, callFqn: String) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowStdoutStderr(id: WorkflowId) extends ApiHandlerMessage
  final case class ApiHandlerCallCaching(id: WorkflowId, parameters: QueryParameters, callName: Option[String]) extends ApiHandlerMessage

  sealed trait WorkflowManagerResponse

  sealed trait WorkflowManagerSuccessResponse extends WorkflowManagerResponse

  sealed trait WorkflowManagerFailureResponse extends WorkflowManagerResponse {
    def failure: Throwable
  }

  final case class WorkflowManagerSubmitSuccess(id: WorkflowId) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerSubmitFailure(override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerWorkflowOutputsSuccess(id: WorkflowId, outputs: core.WorkflowOutputs) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerWorkflowOutputsFailure(id: WorkflowId, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerStatusSuccess(id: WorkflowId, state: WorkflowState) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerStatusFailure(id: WorkflowId, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerAbortSuccess(id: WorkflowId) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerAbortFailure(id: WorkflowId, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerQuerySuccess(uri: Uri, response: WorkflowQueryResponse, meta: Option[QueryMetadata]) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerQueryFailure(override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerCallOutputsSuccess(id: WorkflowId, callFqn: FullyQualifiedName, outputs: core.JobOutputs) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerCallOutputsFailure(id: WorkflowId, callFqn: FullyQualifiedName, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerCallStdoutStderrFailure(id: WorkflowId, callFqn: FullyQualifiedName, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerWorkflowStdoutStderrFailure(id: WorkflowId, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerWorkflowMetadataFailure(id: WorkflowId, override val failure: Throwable) extends WorkflowManagerFailureResponse
  final case class WorkflowManagerCallCachingSuccess(id: WorkflowId, updateCount: Int) extends WorkflowManagerSuccessResponse
  final case class WorkflowManagerCallCachingFailure(id: WorkflowId, override val failure: Throwable) extends WorkflowManagerFailureResponse
}

class CromwellApiHandler(requestHandlerActor: ActorRef) extends Actor with WorkflowQueryPagination {
  import CromwellApiHandler._
  import WorkflowJsonSupport._

  implicit val timeout = Timeout(2 seconds)
  val log = Logging(context.system, classOf[CromwellApiHandler])
  val conf = ConfigFactory.load()

  def callNotFound(callFqn: String, id: WorkflowId) = {
    RequestComplete(StatusCodes.NotFound, APIResponse.error(new Throwable(s"Call $callFqn not found for workflow '$id'.")))
  }

  private def error(t: Throwable)(f: Throwable => RequestComplete[_]): Unit = context.parent ! f(t)

  override def receive = {

    case ApiHandlerWorkflowAbort(id) => requestHandlerActor ! WorkflowManagerActor.AbortWorkflowCommand(id)
    case WorkflowManagerAbortSuccess(id) =>
      context.parent ! RequestComplete(StatusCodes.OK, WorkflowAbortResponse(id.toString, WorkflowAborted.toString))
    case WorkflowManagerAbortFailure(_, e) =>
      error(e) {
        case _: IllegalStateException => RequestComplete(StatusCodes.Forbidden, APIResponse.error(e))
        case _: WorkflowNotFoundException => RequestComplete(StatusCodes.NotFound, APIResponse.error(e))
        case _ => RequestComplete(StatusCodes.InternalServerError, APIResponse.error(e))
      }

    case ApiHandlerWorkflowSubmit(source) => requestHandlerActor ! WorkflowStoreActor.SubmitWorkflow(source)

    case WorkflowStoreActor.WorkflowSubmittedToStore(id) =>
      context.parent ! RequestComplete(StatusCodes.Created, WorkflowSubmitResponse(id.toString, WorkflowSubmitted.toString))

    case ApiHandlerWorkflowSubmitBatch(sources) => requestHandlerActor ! WorkflowStoreActor.BatchSubmitWorkflows(sources)

    case WorkflowStoreActor.WorkflowsBatchSubmittedToStore(ids) =>
      val responses = ids map { id => WorkflowSubmitResponse(id.toString, WorkflowSubmitted.toString) }
      context.parent ! RequestComplete(StatusCodes.OK, responses.list)
  }
}
