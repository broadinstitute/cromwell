package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import akka.actor.ActorRef
import akka.pattern.ask
import akka.util.Timeout
import cats.data.NonEmptyVector
import cromwell.core.WorkflowId

import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{ExecutionContext, Future}

/**
  * Interface for workflow store operations that read or write multiple rows in a single transaction.
  */
sealed trait WorkflowStoreAccess {
  def writeWorkflowHeartbeats(workflowIds: NonEmptyVector[(WorkflowId, OffsetDateTime)],
                              heartbeatDateTime: OffsetDateTime)
                             (implicit ec: ExecutionContext): Future[Int]

  def fetchStartableWorkflows(maxWorkflows: Int, cromwellId: String, heartbeatTtl: FiniteDuration)
                             (implicit ec: ExecutionContext): Future[List[WorkflowToStart]]
}

/**
  * An implementation of `WorkflowStoreAccess` that makes no attempt to coordinate requests to the workflow store
  * and simply forwards them to the underlying database layer.
  */
case class UncoordinatedWorkflowStoreAccess(store: WorkflowStore) extends WorkflowStoreAccess {

  override def writeWorkflowHeartbeats(workflowIds: NonEmptyVector[(WorkflowId, OffsetDateTime)],
                                       heartbeatDateTime: OffsetDateTime)
                                      (implicit ec: ExecutionContext): Future[Int] = {
    store.writeWorkflowHeartbeats(workflowIds.toVector.toSet, heartbeatDateTime)
  }

  override def fetchStartableWorkflows(maxWorkflows: Int, cromwellId: String, heartbeatTtl: FiniteDuration)
                                      (implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
    store.fetchStartableWorkflows(maxWorkflows, cromwellId, heartbeatTtl)
  }
}

/**
  * An implementation of `WorkflowStoreAccess` that coordinates access to the workflow store behind an actor
  * that runs its operations sequentially.
  */
case class CoordinatedWorkflowStoreAccess(actor: ActorRef) extends WorkflowStoreAccess {
  override def writeWorkflowHeartbeats(workflowIds: NonEmptyVector[(WorkflowId, OffsetDateTime)],
                                       heartbeatDateTime: OffsetDateTime)
                                      (implicit ec: ExecutionContext): Future[Int] = {
    implicit val timeout = Timeout(WorkflowStoreCoordinatedAccessActor.Timeout)
    actor.ask(WorkflowStoreCoordinatedAccessActor.WriteHeartbeats(workflowIds, heartbeatDateTime)).mapTo[Int]
  }

  override def fetchStartableWorkflows(maxWorkflows: Int, cromwellId: String, heartbeatTtl: FiniteDuration)
                                      (implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
    implicit val timeout = Timeout(WorkflowStoreCoordinatedAccessActor.Timeout)
    val message = WorkflowStoreCoordinatedAccessActor.FetchStartableWorkflows(maxWorkflows, cromwellId, heartbeatTtl)
    actor.ask(message).mapTo[List[WorkflowToStart]]
  }
}

