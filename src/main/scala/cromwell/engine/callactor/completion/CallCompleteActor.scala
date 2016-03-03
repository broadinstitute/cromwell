package cromwell.engine.callactor.completion

import java.sql.SQLTimeoutException

import akka.actor.{ActorRef, FSM, PoisonPill, Props}
import cromwell.engine.CromwellFSM.ActorFailure
import cromwell.engine._
import cromwell.engine.backend.BackendCall
import cromwell.engine.callactor.completion.CallCompleteActor._
import cromwell.engine.db.DBActor.{PersistCallCompleteData, PersistCallStatusTerminal}
import cromwell.engine.db.DBQueryActor.{FailedQuery, SuccessfulQuery}
import cromwell.engine.db.{DBActor, ExecutionDatabaseKey}
import cromwell.engine.workflow.{ExecutionStoreKey, OutputKey, WorkflowActor}
import cromwell.logging.WorkflowLogger

object CallCompleteActor {
  sealed trait CallCompleteState
  case object Idle extends CallCompleteState
  case object PersistingOutputs extends CallCompleteState {
    def failureBuilder(key: OutputKey)(failure: Throwable) = {
      s"Outputs / Events persistence failed for call ${key.tag}."
    }
  }
  case object PersistingStatus extends CallCompleteState
  case object Done extends CallCompleteState

  sealed trait CallCompleteMessage
  private case object Start extends CallCompleteMessage

  implicit class EnhancedExecutionStoreKey(val key: ExecutionStoreKey) extends AnyVal {
    def toDatabaseKey: ExecutionDatabaseKey = ExecutionDatabaseKey(key.scope.fullyQualifiedName, key.index, key.attempt)
  }

  def props(workflowId: WorkflowId,
            key: OutputKey,
            outputs: CallOutputs,
            events: Seq[ExecutionEventEntry],
            returnCode: Int,
            hash: Option[ExecutionHash],
            clonedFrom: Option[BackendCall],
            replyTo: ActorRef,
            logger: WorkflowLogger) = {
    Props(new CallCompleteActor(workflowId, key, outputs, events, returnCode, hash, clonedFrom, replyTo, logger, DBActor.instance))
  }

}

class CallCompleteActor(workflowId: WorkflowId,
                        key: OutputKey,
                        outputs: CallOutputs,
                        events: Seq[ExecutionEventEntry],
                        returnCode: Int,
                        hash: Option[ExecutionHash],
                        clonedFrom: Option[BackendCall],
                        replyTo: ActorRef,
                        val logger: WorkflowLogger,
                        dbActor: ActorRef) extends CromwellFSM[CallCompleteState, Unit] with RetryableFSM[CallCompleteState, Unit] {

  startWith(Idle, ())
  self ! Start

  when(Idle) {
    case Event(m @ Start, _) =>
      dbActor ! DBActor.PersistCallCompleteData(workflowId, key, outputs, events, logger)
      goto(PersistingOutputs)
  }

  when(PersistingOutputs) {
    case Event(SuccessfulQuery(m: PersistCallCompleteData, _), _) =>
      logger.info(s"persisting status of ${key.tag} to Done.")
      dbActor ! DBActor.PersistCallStatusTerminal(workflowId, key.toDatabaseKey, ExecutionStatus.Done, Option(returnCode), hash, clonedFrom, logger)
      goto(PersistingStatus)
  }

  when(PersistingStatus) {
    case Event(SuccessfulQuery(m: PersistCallStatusTerminal, _), _) =>
      replyTo ! WorkflowActor.CallSucceeded(key, outputs)
      self ! PoisonPill
      goto(Done)
  }

  when(Done) { FSM.NullFunction }

  whenUnhandled {
    case Event(f @ FailedQuery(_, _, message), _)  =>
      retryOrFail(f, message, _ => WorkflowActor.CallCompletionFailed(key))
  }

  override val failureState: CallCompleteState = Done
  override def sendFailureTo: Option[ActorRef] = Option(replyTo)

  override def isRetryable(failure: ActorFailure): Boolean = failure match {
    case _: SQLTimeoutException => true
    case _ => false
  }
}
