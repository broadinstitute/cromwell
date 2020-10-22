package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import akka.actor.ActorRef
import akka.pattern.ask
import akka.util.Timeout
import cats.data.NonEmptyVector
import cromwell.core.WorkflowId
import cromwell.engine.workflow.workflowstore.SqlWorkflowStore.WorkflowStoreAbortResponse.WorkflowStoreAbortResponse

import scala.concurrent.duration.FiniteDuration
import scala.concurrent.{ExecutionContext, Future}

/**
  * Interface for workflow store operations that empirically come into conflict and deadlock [WA-334]
  *
  * The problematic pairs all involve fetching:
  * fetch <-> heartbeat
  * fetch <-> abort
  * fetch <-> delete
  *
  */
sealed trait WorkflowStoreAccess {
  def writeWorkflowHeartbeats(workflowIds: NonEmptyVector[(WorkflowId, OffsetDateTime)],
                              heartbeatDateTime: OffsetDateTime)
                             (implicit ec: ExecutionContext): Future[Int]

  def fetchStartableWorkflows(maxWorkflows: Int, cromwellId: String, heartbeatTtl: FiniteDuration)
                             (implicit ec: ExecutionContext): Future[List[WorkflowToStart]]

  def abort(workflowId: WorkflowId)
           (implicit ec: ExecutionContext): Future[WorkflowStoreAbortResponse]

  def deleteFromStore(workflowId: WorkflowId)
                     (implicit ec: ExecutionContext): Future[Int]

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

  override def deleteFromStore(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Int] = {
    store.deleteFromStore(workflowId)
  }

  override def abort(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowStoreAbortResponse] = {
    store.abort(workflowId)
  }
}

/**
  * An implementation of `WorkflowStoreAccess` that coordinates access to the workflow store behind an actor
  * that runs its operations sequentially. Enabled by default in `CromwellRootActor`.
  */
case class CoordinatedWorkflowStoreAccess(coordinatedWorkflowStoreAccessActor: ActorRef) extends WorkflowStoreAccess {
  implicit val timeout = Timeout(WorkflowStoreCoordinatedAccessActor.Timeout)

  override def writeWorkflowHeartbeats(workflowIds: NonEmptyVector[(WorkflowId, OffsetDateTime)],
                                       heartbeatDateTime: OffsetDateTime)
                                      (implicit ec: ExecutionContext): Future[Int] = {
    coordinatedWorkflowStoreAccessActor.ask(WorkflowStoreCoordinatedAccessActor.WriteHeartbeats(workflowIds, heartbeatDateTime)).mapTo[Int]
  }

  override def fetchStartableWorkflows(maxWorkflows: Int, cromwellId: String, heartbeatTtl: FiniteDuration)
                                      (implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = {
    val message = WorkflowStoreCoordinatedAccessActor.FetchStartableWorkflows(maxWorkflows, cromwellId, heartbeatTtl)
    coordinatedWorkflowStoreAccessActor.ask(message).mapTo[List[WorkflowToStart]]
  }

  override def deleteFromStore(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Int] = {
    coordinatedWorkflowStoreAccessActor.ask(WorkflowStoreCoordinatedAccessActor.DeleteFromStore(workflowId)).mapTo[Int]
  }

  override def abort(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowStoreAbortResponse] = {
    coordinatedWorkflowStoreAccessActor.ask(WorkflowStoreCoordinatedAccessActor.Abort(workflowId)).mapTo[WorkflowStoreAbortResponse]
  }
}

