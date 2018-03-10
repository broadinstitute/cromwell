package cromwell.engine.workflow.workflowstore

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cats.data.NonEmptyList
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core._
import cromwell.util.GracefulShutdownHelper
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.concurrent.duration.FiniteDuration

final case class WorkflowStoreActor private(store: WorkflowStore, serviceRegistryActor: ActorRef, abortAllJobsOnTerminate: Boolean, cromwellId: String, heartbeatTtl: FiniteDuration)
  extends Actor with ActorLogging with GracefulShutdownHelper {
  import WorkflowStoreActor._

  lazy val workflowStoreSubmitActor: ActorRef = context.actorOf(WorkflowStoreSubmitActor.props(store, serviceRegistryActor), "WorkflowStoreSubmitActor")
  lazy val workflowStoreEngineActor: ActorRef = context.actorOf(
    WorkflowStoreEngineActor.props(store, serviceRegistryActor, abortAllJobsOnTerminate, cromwellId, heartbeatTtl), "WorkflowStoreEngineActor")

  override def receive = {
    case ShutdownCommand => waitForActorsAndShutdown(NonEmptyList.of(workflowStoreSubmitActor, workflowStoreEngineActor))
    case cmd: WorkflowStoreActorSubmitCommand => workflowStoreSubmitActor forward cmd
    case cmd: WorkflowStoreActorEngineCommand => workflowStoreEngineActor forward cmd
  }
}

object WorkflowStoreActor {
  sealed trait WorkflowStoreActorEngineCommand
  final case class FetchRunnableWorkflows(n: Int) extends WorkflowStoreActorEngineCommand
  final case class AbortWorkflowCommand(id: WorkflowId) extends WorkflowStoreActorEngineCommand
  case object InitializerCommand extends WorkflowStoreActorEngineCommand
  case object WorkDone extends WorkflowStoreActorEngineCommand
  case object AbortAllRunningWorkflowsCommandAndStop extends WorkflowStoreActorEngineCommand

  sealed trait WorkflowStoreActorSubmitCommand
  final case class SubmitWorkflow(source: WorkflowSourceFilesCollection) extends WorkflowStoreActorSubmitCommand
  final case class BatchSubmitWorkflows(sources: NonEmptyList[WorkflowSourceFilesCollection]) extends WorkflowStoreActorSubmitCommand

  def props(workflowStoreDatabase: WorkflowStore, serviceRegistryActor: ActorRef, abortAllJobsOnTerminate: Boolean, cromwellId: String, heartbeatTtl: FiniteDuration) = {
    Props(WorkflowStoreActor(workflowStoreDatabase, serviceRegistryActor, abortAllJobsOnTerminate, cromwellId, heartbeatTtl)).withDispatcher(EngineDispatcher)
  }
}
