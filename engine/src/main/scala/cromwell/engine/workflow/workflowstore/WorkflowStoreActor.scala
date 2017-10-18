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
    case WorkflowStatusRequest(id) => statusOf(id, sender())
    case ValidateWorkflowId(id) => checkWorkflowIdExistence(id, sender())
    case cmd: WorkflowStoreActorSubmitCommand => workflowStoreSubmitActor forward cmd
    case cmd: WorkflowStoreActorEngineCommand => workflowStoreEngineActor forward cmd
  }

  private def statusOf(id: WorkflowId, replyTo: ActorRef): Unit = {
    pipe(store.status(id)).to(replyTo)
    ()
  }

  private def checkWorkflowIdExistence(id: WorkflowId, replyTo: ActorRef): Unit = {
    val validation: Future[WorkflowValidationResponse] = store.status(id) map {
      case Some(_) => RecognizedWorkflowId
      case None => UnrecognizedWorkflowId
    } recover {
      case failure => FailedToCheckWorkflowId(failure)
    }

    pipe(validation).to(replyTo)
    ()
  }
}

object WorkflowStoreActor {
  sealed trait WorkflowStoreActorEngineCommand
  final case class FetchRunnableWorkflows(n: Int) extends WorkflowStoreActorEngineCommand
  final case class AbortWorkflow(id: WorkflowId, manager: ActorRef) extends WorkflowStoreActorEngineCommand
  case object InitializerCommand extends WorkflowStoreActorEngineCommand
  case object WorkDone extends WorkflowStoreActorEngineCommand

  sealed trait WorkflowStoreActorSubmitCommand
  final case class SubmitWorkflow(source: WorkflowSourceFilesCollection) extends WorkflowStoreActorSubmitCommand
  final case class BatchSubmitWorkflows(sources: NonEmptyList[WorkflowSourceFilesCollection]) extends WorkflowStoreActorSubmitCommand
  final case class WorkflowStatusRequest(id: WorkflowId) extends WorkflowStoreActorSubmitCommand

  def props(workflowStoreDatabase: WorkflowStore, serviceRegistryActor: ActorRef) = {
    Props(WorkflowStoreActor(workflowStoreDatabase, serviceRegistryActor)).withDispatcher(EngineDispatcher)
  }
}
