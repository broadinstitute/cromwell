package cromwell.webservice

import akka.actor.{Actor, ActorRef, Props}
import akka.event.Logging
import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.core._
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.engine.workflow.WorkflowManagerActor.WorkflowNotFoundException
import cromwell.engine.workflow.workflowstore.{WorkflowStoreActor, WorkflowStoreEngineActor, WorkflowStoreSubmitActor}
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.webservice.metadata.WorkflowQueryPagination
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._


object CromwellApiHandler {
  def props(requestHandlerActor: ActorRef): Props = {
    Props(new CromwellApiHandler(requestHandlerActor)).withDispatcher(ApiDispatcher)
  }

  sealed trait ApiHandlerMessage

  final case class ApiHandlerWorkflowSubmit(source: WorkflowSourceFilesCollection) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowSubmitBatch(sources: NonEmptyList[WorkflowSourceFilesCollection]) extends ApiHandlerMessage
  final case class ApiHandlerWorkflowAbort(id: WorkflowId, manager: ActorRef) extends ApiHandlerMessage
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
    case WorkflowStoreEngineActor.WorkflowAborted(id) =>
      context.parent ! RequestComplete((StatusCodes.OK, WorkflowAbortResponse(id.toString, WorkflowAborted.toString)))
    case WorkflowStoreEngineActor.WorkflowAbortFailed(_, e) =>
      error(e) {
        case _: IllegalStateException => RequestComplete((StatusCodes.Forbidden, APIResponse.error(e)))
        case _: WorkflowNotFoundException => RequestComplete((StatusCodes.NotFound, APIResponse.error(e)))
        case _ => RequestComplete((StatusCodes.InternalServerError, APIResponse.error(e)))
      }

    case ApiHandlerWorkflowSubmit(source) => requestHandlerActor ! WorkflowStoreActor.SubmitWorkflow(source)

    case WorkflowStoreSubmitActor.WorkflowSubmittedToStore(id) =>
      context.parent ! RequestComplete((StatusCodes.Created, WorkflowSubmitResponse(id.toString, WorkflowSubmitted.toString)))

    case ApiHandlerWorkflowSubmitBatch(sources) => requestHandlerActor !
      WorkflowStoreActor.BatchSubmitWorkflows(sources.map(w =>
        WorkflowSourceFilesCollection(
          workflowSource = w.workflowSource,
          workflowType = w.workflowType,
          workflowTypeVersion = w.workflowTypeVersion,
          inputsJson = w.inputsJson,
          workflowOptionsJson = w.workflowOptionsJson,
          labelsJson = w.labelsJson,
          importsFile = w.importsZipFileOption)))


    case WorkflowStoreSubmitActor.WorkflowsBatchSubmittedToStore(ids) =>
      val responses = ids map { id => WorkflowSubmitResponse(id.toString, WorkflowSubmitted.toString) }
      context.parent ! RequestComplete((StatusCodes.OK, responses.toList))
  }
}
