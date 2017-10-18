package cromwell.engine.workflow.workflowstore

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.pattern.pipe
import cats.data.NonEmptyList
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core._
import cromwell.services.metadata.MetadataService._
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.Future

final case class WorkflowStoreActor private(store: WorkflowStore, serviceRegistryActor: ActorRef)
  extends Actor with ActorLogging with GracefulShutdownHelper {
  import WorkflowStoreActor._
  implicit val ec = context.dispatcher
  
  lazy val workflowStoreSubmitActor: ActorRef = context.actorOf(WorkflowStoreSubmitActor.props(store, serviceRegistryActor), "WorkflowStoreSubmitActor")
  lazy val workflowStoreEngineActor: ActorRef = context.actorOf(
    WorkflowStoreEngineActor.props(store, serviceRegistryActor), "WorkflowStoreEngineActor")

  override def receive = {
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(workflowStoreSubmitActor))
    case WorkflowStateRequest(id) => 
      pipe(statusOf(id)).to(sender())
      ()
    case ValidateWorkflowId(id) => 
      pipe(checkWorkflowIdExistence(id)).to(sender())
      ()
    case cmd: WorkflowStoreActorSubmitCommand => workflowStoreSubmitActor forward cmd
    case cmd: WorkflowStoreActorEngineCommand => workflowStoreEngineActor forward cmd
  }

  private def statusOf(id: WorkflowId): Future[WorkflowStateResponse] = {
    store.status(id) map {
      case Some(state) => WorkflowStateSuccessfulResponse(id, state)
      case None => WorkflowStateNotFoundResponse(id)
    } recover {
      case failure => WorkflowStateFailedResponse(id, failure)
    }
  }

  private def checkWorkflowIdExistence(id: WorkflowId): Future[WorkflowValidationResponse] = {
    statusOf(id).map(_.toValidationResponse)
  }
}

object WorkflowStoreActor {
  sealed trait WorkflowStoreCommand
  sealed trait WorkflowStoreActorEngineCommand extends WorkflowStoreCommand
  final case class FetchRunnableWorkflows(n: Int) extends WorkflowStoreActorEngineCommand
  final case class AbortWorkflow(id: WorkflowId, manager: ActorRef) extends WorkflowStoreActorEngineCommand
  case object InitializerCommand extends WorkflowStoreActorEngineCommand
  case object WorkDone extends WorkflowStoreActorEngineCommand

  sealed trait WorkflowStoreActorSubmitCommand extends WorkflowStoreCommand
  final case class SubmitWorkflow(source: WorkflowSourceFilesCollection) extends WorkflowStoreActorSubmitCommand
  final case class BatchSubmitWorkflows(sources: NonEmptyList[WorkflowSourceFilesCollection]) extends WorkflowStoreActorSubmitCommand
  
  final case class WorkflowStateRequest(id: WorkflowId) extends WorkflowStoreCommand

  sealed trait WorkflowStateResponse {
    def toValidationResponse: WorkflowValidationResponse
  }
  final case class WorkflowStateSuccessfulResponse(id: WorkflowId, status: WorkflowStoreState) extends WorkflowStateResponse {
    override def toValidationResponse: WorkflowValidationResponse = RecognizedWorkflowId
  }
  final case class WorkflowStateNotFoundResponse(id: WorkflowId) extends WorkflowStateResponse {
    override def toValidationResponse: WorkflowValidationResponse = UnrecognizedWorkflowId
  }
  final case class WorkflowStateFailedResponse(id: WorkflowId, ex: Throwable) extends WorkflowStateResponse {
    override def toValidationResponse: WorkflowValidationResponse = FailedToCheckWorkflowId(ex)
  }

  def props(workflowStoreDatabase: WorkflowStore, serviceRegistryActor: ActorRef) = {
    Props(WorkflowStoreActor(workflowStoreDatabase, serviceRegistryActor)).withDispatcher(EngineDispatcher)
  }
}
