package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import akka.actor.{Actor, ActorSystem, Props, Status}
import cats.data.NonEmptyVector
import cromwell.core.{Dispatcher, WorkflowId}
import cromwell.engine.workflow.workflowstore.WorkflowStoreCoordinatedAccessActor._
import mouse.all._

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

/**
  * Serializes access to the workflow store for workflow store writers that acquire locks to multiple rows inside a single
  * transaction and otherwise are prone to deadlock.
  */
class WorkflowStoreCoordinatedAccessActor(workflowStore: WorkflowStore) extends Actor {
  implicit val ec: ExecutionContext = context.system.dispatcher
  implicit val actorSystem: ActorSystem = context.system

  def run[A](future: Future[A]): Unit = {
    val result = Try(Await.result(future, Timeout)) match {
      case Success(s) => s
      case f: Failure[_] => Status.Failure(f.exception)
    }
    sender() ! result
  }

  override def receive: Receive = {
    case WriteHeartbeats(ids, heartbeatDateTime) =>
      workflowStore.writeWorkflowHeartbeats(ids.toVector.toSet, heartbeatDateTime) |> run
    case FetchStartableWorkflows(count, cromwellId, heartbeatTtl, excludedGroups) =>
      workflowStore.fetchStartableWorkflows(count, cromwellId, heartbeatTtl, excludedGroups) |> run
    case DeleteFromStore(workflowId) =>
      workflowStore.deleteFromStore(workflowId) |> run
    case Abort(workflowId) =>
      workflowStore.abort(workflowId) |> run
  }
}

object WorkflowStoreCoordinatedAccessActor {
  final case class WriteHeartbeats(workflowIds: NonEmptyVector[(WorkflowId, OffsetDateTime)],
                                   heartbeatDateTime: OffsetDateTime
  )
  final case class FetchStartableWorkflows(count: Int,
                                           cromwellId: String,
                                           heartbeatTtl: FiniteDuration,
                                           excludedGroups: Set[String]
  )
  final case class DeleteFromStore(workflowId: WorkflowId)
  final case class Abort(workflowId: WorkflowId)

  val Timeout = 1 minute

  def props(workflowStore: WorkflowStore): Props =
    Props(new WorkflowStoreCoordinatedAccessActor(workflowStore)).withDispatcher(Dispatcher.IoDispatcher)
}
