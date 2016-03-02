package cromwell.engine.callactor.completion

import akka.actor.{FSM, ActorRef, PoisonPill, Props}
import cromwell.engine._
import cromwell.engine.backend.BackendCall
import cromwell.engine.callactor.completion.CallCompleteActor._
import cromwell.engine.db.DataAccess._
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.workflow.{WorkflowActor, ExecutionStoreKey, OutputKey}
import cromwell.logging.WorkflowLogger

import scala.concurrent.ExecutionContext.Implicits.global

object CallCompleteActor {
  sealed trait CallCompleteState
  case object Idle extends CallCompleteState
  case object PersistingOutputs extends CallCompleteState
  case object PersistingStatus extends CallCompleteState
  case object Done extends CallCompleteState

  sealed trait CallCompleteMessage
  private case object Start extends CallCompleteMessage
  private case object OutputsPersisted extends CallCompleteMessage
  private case object StatusPersisted extends CallCompleteMessage

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
    Props(new CallCompleteActor(workflowId, key, outputs, events, returnCode, hash, clonedFrom, replyTo, logger))
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
                        val logger: WorkflowLogger) extends CromwellFSM[CallCompleteState, Unit] {

  startWith(Idle, ())
  self ! Start

  when(Idle) {
    case Event(Start, _) =>
      val completionWork = for {
        _ <- globalDataAccess.setOutputs(workflowId, key, outputs, key.scope.rootWorkflow.outputs)
        _ <- globalDataAccess.setExecutionEvents(workflowId, key.scope.fullyQualifiedName, key.index, key.attempt, events)
      } yield OutputsPersisted

      chain(completionWork, PersistingOutputs, s"Outputs / Events persistence failed for call ${key.tag}.")
  }

  when(PersistingOutputs) {
    case Event(OutputsPersisted, _) =>
      val result = globalDataAccess.setTerminalStatus(workflowId, key.toDatabaseKey, ExecutionStatus.Done, Option(returnCode), hash, clonedFrom) map {
        _ => StatusPersisted
      }

      chain(result, PersistingStatus, s"Stats persistence failed for call ${key.tag}.")
  }

  when(PersistingStatus) {
    case Event(StatusPersisted, _) =>
      replyTo ! WorkflowActor.CallSucceeded(key, outputs)
      self ! PoisonPill
      goto(Done)
  }

  when(Done) { FSM.NullFunction }

  whenUnhandled {
    case Event(f @ ActorFailure(_, _), _)  =>
      fail(f, WorkflowActor.CallCompletionFailed(key))
  }

  override val failureState: CallCompleteState = Done
  override def sendFailureTo: Option[ActorRef] = Option(replyTo)
}
