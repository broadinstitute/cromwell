package cromwell.engine.callactor.completion

import akka.actor.{FSM, ActorRef, PoisonPill, Props}
import cromwell.engine._
import cromwell.engine.callactor.completion.CallRetriedActor._
import cromwell.engine.db.DataAccess._
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.workflow.{ExecutionStoreKey, OutputKey, WorkflowActor}
import cromwell.logging.WorkflowLogger

import scala.concurrent.ExecutionContext.Implicits.global

object CallRetriedActor {
  sealed trait CallRetriedState
  case object Idle extends CallRetriedState
  case object PersistingStatus extends CallRetriedState
  case object Done extends CallRetriedState

  sealed trait CallRetriedMessage
  private case object Start extends CallRetriedMessage
  private case object StatusPersisted extends CallRetriedMessage

  implicit class EnhancedExecutionStoreKey(val key: ExecutionStoreKey) extends AnyVal {
    def toDatabaseKey: ExecutionDatabaseKey = ExecutionDatabaseKey(key.scope.fullyQualifiedName, key.index, key.attempt)
  }

  def props(workflowId: WorkflowId,
            key: OutputKey,
            events: Seq[ExecutionEventEntry],
            returnCode: Option[Int],
            replyTo: ActorRef,
            logger: WorkflowLogger) = {
    Props(new CallRetriedActor(workflowId, key, events, returnCode, replyTo, logger))
  }

}

class CallRetriedActor(workflowId: WorkflowId,
                        key: OutputKey,
                        events: Seq[ExecutionEventEntry],
                        returnCode: Option[Int],
                        replyTo: ActorRef,
                        val logger: WorkflowLogger) extends CromwellFSM[CallRetriedState, Unit] {

  startWith(Idle, ())
  self ! Start

  when(Idle) {
    case Event(Start, _) =>
      val completionWork = for {
        _ <- globalDataAccess.setTerminalStatus(workflowId, key.toDatabaseKey, ExecutionStatus.Preempted, returnCode, None, None)
        _ <- globalDataAccess.setExecutionEvents(workflowId, key.scope.fullyQualifiedName, key.index, key.attempt, events)
      } yield StatusPersisted

      chain(completionWork, PersistingStatus, s"Outputs / Events persistence failed for call ${key.tag}.")
  }

  when(PersistingStatus) {
    case Event(StatusPersisted, _) =>
      replyTo ! WorkflowActor.CallFailedRetryable(key)
      self ! PoisonPill
      goto(Done)
  }

  when(Done) { FSM.NullFunction }

  whenUnhandled {
    case Event(f @ ActorFailure(_, _), _)  =>
      fail(f, WorkflowActor.CallCompletionFailed(key))
  }

  override val failureState: CallRetriedState = Done
  override def sendFailureTo: Option[ActorRef] = Option(replyTo)
}
