package cromwell.webservice

import akka.actor.{Actor, ActorRef, Props}
import akka.event.Logging
import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.core._
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.webservice.metadata.WorkflowQueryPagination
import spray.http.{StatusCodes, Uri}
import spray.httpx.SprayJsonSupport._


object CromwellApiHandler {
  def props(requestHandlerActor: ActorRef): Props = {
    Props(new CromwellApiHandler(requestHandlerActor)).withDispatcher(ApiDispatcher)
  }

  sealed trait ApiHandlerMessage

  final case class ApiHandlerWorkflowSubmit(source: WorkflowSourceFilesCollection) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowSubmitBatch(sources: NonEmptyList[WorkflowSourceFilesCollection]) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowQuery(uri: Uri, parameters: Seq[(String, String)]) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowStatus(id: WorkflowId) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowOutputs(id: WorkflowId) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowAbort(id: WorkflowId, manager: ActorRef) extends ApiHandlerMessage
  final case class ApiHandlerCallOutputs(id: WorkflowId, callFqn: String) extends ApiHandlerMessage
  final case class ApiHandlerCallStdoutStderr(id: WorkflowId, callFqn: String) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowStdoutStderr(id: WorkflowId) extends ApiHandlerMessage
  final case class ApiHandlerCallCaching(id: WorkflowId, parameters: QueryParameters, callName: Option[String]) extends ApiHandlerMessage
  case object ApiHandlerEngineStats extends ApiHandlerMessage
}

class CromwellApiHandler(requestHandlerActor: ActorRef) extends Actor with WorkflowQueryPagination {
  import CromwellApiHandler._
  import WorkflowJsonSupport._

  val log = Logging(context.system, classOf[CromwellApiHandler])
  val conf = ConfigFactory.load()

  def callNotFound(callFqn: String, id: WorkflowId) = {
    RequestComplete((StatusCodes.NotFound, APIResponse.error(
      new RuntimeException(s"Call $callFqn not found for workflow '$id'."))))
  }

  private def error(t: Throwable)(f: Throwable => RequestComplete[_]): Unit = context.parent ! f(t)

  override def receive = {
    case ApiHandlerEngineStats => requestHandlerActor ! WorkflowManagerActor.EngineStatsCommand
    case stats: EngineStatsActor.EngineStats => context.parent ! RequestComplete((StatusCodes.OK, stats))
    case ApiHandlerWorkflowAbort(id, manager) => requestHandlerActor ! WorkflowStoreActor.AbortWorkflow(id, manager)
    case WorkflowStoreActor.WorkflowAborted(id) =>
      context.parent ! RequestComplete((StatusCodes.OK, WorkflowAbortResponse(id.toString, WorkflowAborted.toString)))
    case WorkflowStoreActor.WorkflowAbortFailed(_, e) =>
      error(e) {
        case _: IllegalStateException => RequestComplete((StatusCodes.Forbidden, APIResponse.error(e)))
        case _: WorkflowNotFoundException => RequestComplete((StatusCodes.NotFound, APIResponse.error(e)))
        case _ => RequestComplete((StatusCodes.InternalServerError, APIResponse.error(e)))
      }

    case ApiHandlerWorkflowSubmit(source) => requestHandlerActor ! WorkflowStoreActor.SubmitWorkflow(source)

    case WorkflowStoreActor.WorkflowSubmittedToStore(id) =>
      context.parent ! RequestComplete((StatusCodes.Created, WorkflowSubmitResponse(id.toString, WorkflowSubmitted.toString)))

    case ApiHandlerWorkflowSubmitBatch(sources) => requestHandlerActor !
      WorkflowStoreActor.BatchSubmitWorkflows(sources.map(x => WorkflowSourceFilesCollection(x.wdlSource,x.inputsJson,x.workflowOptionsJson,x.importsZipFileOption)))


    case WorkflowStoreActor.WorkflowsBatchSubmittedToStore(ids) =>
      val responses = ids map { id => WorkflowSubmitResponse(id.toString, WorkflowSubmitted.toString) }
      context.parent ! RequestComplete((StatusCodes.OK, responses.toList))
  }
}
