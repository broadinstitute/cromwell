package cromwell.engine.workflow.workflowstore

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.core._
import cromwell.core.Dispatcher.EngineDispatcher

final case class WorkflowStoreActor(store: WorkflowStore, serviceRegistryActor: ActorRef) extends Actor with ActorLogging {
  import WorkflowStoreActor._

  lazy val workflowStoreSubmitActor: ActorRef = context.actorOf(WorkflowStoreSubmitActor.props(store, serviceRegistryActor))
  lazy val workflowStoreEngineActor: ActorRef = context.actorOf(WorkflowStoreEngineActor.props(store, serviceRegistryActor))

  override def receive = {
    case cmd: WorkflowStoreActorSubmitCommand => workflowStoreSubmitActor forward cmd
    case cmd: WorkflowStoreActorEngineCommand => workflowStoreEngineActor forward cmd
  }
}

object WorkflowStoreActor {
  sealed trait WorkflowStoreActorEngineCommand
  final case class FetchRunnableWorkflows(n: Int) extends WorkflowStoreActorEngineCommand
  final case class RemoveWorkflow(id: WorkflowId) extends WorkflowStoreActorEngineCommand
  final case class AbortWorkflow(id: WorkflowId, manager: ActorRef) extends WorkflowStoreActorEngineCommand
  case object InitializerCommand extends WorkflowStoreActorEngineCommand
  case object WorkDone extends WorkflowStoreActorEngineCommand

  sealed trait WorkflowStoreActorSubmitCommand
  final case class SubmitWorkflow(source: WorkflowSourceFilesCollection) extends WorkflowStoreActorSubmitCommand
  final case class BatchSubmitWorkflows(sources: NonEmptyList[WorkflowSourceFilesCollection]) extends WorkflowStoreActorSubmitCommand

  def props(workflowStoreDatabase: WorkflowStore, serviceRegistryActor: ActorRef) = {
    Props(WorkflowStoreActor(workflowStoreDatabase, serviceRegistryActor)).withDispatcher(EngineDispatcher)
  }
}
