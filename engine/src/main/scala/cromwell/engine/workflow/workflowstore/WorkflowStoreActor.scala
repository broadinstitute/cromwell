package cromwell.engine.workflow.workflowstore

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core._
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

final case class WorkflowStoreActor private(
                                             workflowStore: WorkflowStore,
                                             workflowStoreCoordinatedWriteActor: ActorRef,
                                             serviceRegistryActor: ActorRef,
                                             abortAllJobsOnTerminate: Boolean,
                                             workflowHeartbeatConfig: WorkflowHeartbeatConfig)
  extends Actor with ActorLogging with GracefulShutdownHelper {
  import WorkflowStoreActor._

  lazy val workflowStoreSubmitActor: ActorRef = context.actorOf(
    WorkflowStoreSubmitActor.props(
      workflowStoreDatabase = workflowStore,
      serviceRegistryActor = serviceRegistryActor),
    "WorkflowStoreSubmitActor")

  lazy val workflowStoreEngineActor: ActorRef = context.actorOf(
    WorkflowStoreEngineActor.props(
      workflowStore = workflowStore,
      workflowStoreCoordinatedWriteActor = workflowStoreCoordinatedWriteActor,
      serviceRegistryActor = serviceRegistryActor,
      abortAllJobsOnTerminate = abortAllJobsOnTerminate,
      workflowHeartbeatConfig = workflowHeartbeatConfig),
    "WorkflowStoreEngineActor")

  lazy val workflowStoreHeartbeatWriteActor: ActorRef = context.actorOf(
    WorkflowStoreHeartbeatWriteActor.props(
      workflowStoreCoordindatedWriteActor = workflowStoreCoordinatedWriteActor,
      workflowHeartbeatConfig = workflowHeartbeatConfig,
      serviceRegistryActor = serviceRegistryActor),
    "WorkflowStoreHeartbeatWriteActor")

  override def receive = {
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(workflowStoreSubmitActor, workflowStoreEngineActor))
    case cmd: WorkflowStoreActorSubmitCommand => workflowStoreSubmitActor forward cmd
    case cmd: WorkflowStoreActorEngineCommand => workflowStoreEngineActor forward cmd
    case cmd: WorkflowStoreWriteHeartbeatCommand => workflowStoreHeartbeatWriteActor forward cmd
  }
}

object WorkflowStoreActor {
  sealed trait WorkflowStoreActorEngineCommand
  final case class FetchRunnableWorkflows(n: Int) extends WorkflowStoreActorEngineCommand
  final case class AbortWorkflowCommand(id: WorkflowId) extends WorkflowStoreActorEngineCommand
  final case class WorkflowOnHoldToSubmittedCommand(id: WorkflowId) extends WorkflowStoreActorEngineCommand
  case object InitializerCommand extends WorkflowStoreActorEngineCommand
  case object WorkDone extends WorkflowStoreActorEngineCommand
  case object AbortAllRunningWorkflowsCommandAndStop extends WorkflowStoreActorEngineCommand

  sealed trait WorkflowStoreActorSubmitCommand
  final case class SubmitWorkflow(source: WorkflowSourceFilesCollection) extends WorkflowStoreActorSubmitCommand
  final case class BatchSubmitWorkflows(sources: NonEmptyList[WorkflowSourceFilesCollection]) extends WorkflowStoreActorSubmitCommand

  case class WorkflowStoreWriteHeartbeatCommand(workflowId: WorkflowId)

  def props(
             workflowStoreDatabase: WorkflowStore,
             workflowStoreCoordinatedWriteActor: ActorRef,
             serviceRegistryActor: ActorRef,
             abortAllJobsOnTerminate: Boolean,
             workflowHeartbeatConfig: WorkflowHeartbeatConfig
      ) = {
    Props(WorkflowStoreActor(
      workflowStore = workflowStoreDatabase,
      workflowStoreCoordinatedWriteActor = workflowStoreCoordinatedWriteActor,
      serviceRegistryActor = serviceRegistryActor,
      abortAllJobsOnTerminate = abortAllJobsOnTerminate,
      workflowHeartbeatConfig = workflowHeartbeatConfig)).withDispatcher(EngineDispatcher)
  }
}
