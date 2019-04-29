package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.pattern.pipe
import cats.data.NonEmptyList
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core._
import cromwell.engine.CromwellTerminator
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

final case class WorkflowStoreActor private(
                                             workflowStore: WorkflowStore,
                                             workflowStoreAccess: WorkflowStoreAccess,
                                             serviceRegistryActor: ActorRef,
                                             terminator: CromwellTerminator,
                                             abortAllJobsOnTerminate: Boolean,
                                             workflowHeartbeatConfig: WorkflowHeartbeatConfig)
  extends Actor with ActorLogging with GracefulShutdownHelper {
  import WorkflowStoreActor._

  implicit val ec = context.dispatcher

  lazy val workflowStoreSubmitActor: ActorRef = context.actorOf(
    WorkflowStoreSubmitActor.props(
      workflowStoreDatabase = workflowStore,
      serviceRegistryActor = serviceRegistryActor),
    "WorkflowStoreSubmitActor")

  lazy val workflowStoreEngineActor: ActorRef = context.actorOf(
    WorkflowStoreEngineActor.props(
      workflowStore = workflowStore,
      workflowStoreAccess = workflowStoreAccess,
      serviceRegistryActor = serviceRegistryActor,
      abortAllJobsOnTerminate = abortAllJobsOnTerminate,
      workflowHeartbeatConfig = workflowHeartbeatConfig),
    "WorkflowStoreEngineActor")

  lazy val workflowStoreHeartbeatWriteActor: ActorRef = context.actorOf(
    WorkflowStoreHeartbeatWriteActor.props(
      workflowStoreAccess = workflowStoreAccess,
      workflowHeartbeatConfig = workflowHeartbeatConfig,
      terminator = terminator,
      serviceRegistryActor = serviceRegistryActor),
    "WorkflowStoreHeartbeatWriteActor")

  override def receive = {
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(workflowStoreSubmitActor, workflowStoreEngineActor))
    case cmd: WorkflowStoreActorSubmitCommand => workflowStoreSubmitActor forward cmd
    case cmd: WorkflowStoreActorEngineCommand => workflowStoreEngineActor forward cmd
    case cmd: WorkflowStoreWriteHeartbeatCommand => workflowStoreHeartbeatWriteActor forward cmd
    case GetWorkflowStoreStats =>
      // Retrieve the workflow store stats, convert the WorkflowStoreStates to WorkflowStates
      val stats = workflowStore.stats.map(m => m.map(e => WorkflowState.withName(e._1.toString) -> e._2))
      stats pipeTo sender
      ()
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
  final case class FindWorkflowsWithAbortRequested(cromwellId: String) extends WorkflowStoreActorEngineCommand

  sealed trait WorkflowStoreActorEngineResponse
  case class FindWorkflowsWithAbortRequestedSuccess(ids: Iterable[WorkflowId]) extends WorkflowStoreActorEngineResponse
  case class FindWorkflowsWithAbortRequestedFailure(t: Throwable) extends WorkflowStoreActorEngineResponse

  sealed trait WorkflowStoreActorSubmitCommand
  final case class SubmitWorkflow(source: WorkflowSourceFilesCollection) extends WorkflowStoreActorSubmitCommand
  final case class BatchSubmitWorkflows(sources: NonEmptyList[WorkflowSourceFilesCollection]) extends WorkflowStoreActorSubmitCommand

  final case object GetWorkflowStoreStats

  case class WorkflowStoreWriteHeartbeatCommand(workflowId: WorkflowId, submissionTime: OffsetDateTime)

  def props(
             workflowStoreDatabase: WorkflowStore,
             workflowStoreAccess: WorkflowStoreAccess,
             serviceRegistryActor: ActorRef,
             terminator: CromwellTerminator,
             abortAllJobsOnTerminate: Boolean,
             workflowHeartbeatConfig: WorkflowHeartbeatConfig
      ) = {
    Props(WorkflowStoreActor(
      workflowStore = workflowStoreDatabase,
      workflowStoreAccess = workflowStoreAccess,
      serviceRegistryActor = serviceRegistryActor,
      terminator = terminator,
      abortAllJobsOnTerminate = abortAllJobsOnTerminate,
      workflowHeartbeatConfig = workflowHeartbeatConfig)).withDispatcher(EngineDispatcher)
  }
}
