package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import akka.actor.{Actor, Props, Status}
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
    case FetchStartableWorkflows(count, cromwellId, heartbeatTtl) =>
      workflowStore.fetchStartableWorkflows(count, cromwellId, heartbeatTtl) |> run
  }
}

object WorkflowStoreCoordinatedAccessActor {
  final case class WriteHeartbeats(workflowIds: NonEmptyVector[(WorkflowId, OffsetDateTime)],
                                   heartbeatDateTime: OffsetDateTime)
  final case class FetchStartableWorkflows(count: Int, cromwellId: String, heartbeatTtl: FiniteDuration)

  val Timeout = 1 minute

  def props(workflowStore: WorkflowStore): Props = Props(new WorkflowStoreCoordinatedAccessActor(workflowStore)).withDispatcher(Dispatcher.IoDispatcher)
}
