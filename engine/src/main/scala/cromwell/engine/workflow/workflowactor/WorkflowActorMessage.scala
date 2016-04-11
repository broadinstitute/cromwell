package cromwell.engine.workflow.workflowactor

import akka.actor.ActorRef
import cromwell.core._
import cromwell.engine.ExecutionStatus._
import cromwell.engine.callactor.CallActor
import cromwell.engine.{WorkflowState, ExecutionEventEntry}
import cromwell.engine.backend.{BackendCallJobDescriptor, ExecutionHash, CallLogs}
import cromwell.engine.workflow._

object WorkflowActorMessages {

  sealed trait WorkflowActorMessage

  case object GetFailureMessage extends WorkflowActorMessage

  case object AbortWorkflow extends WorkflowActorMessage

  sealed trait CallMessage extends WorkflowActorMessage {
    def callKey: ExecutionStoreKey
  }

  case class CallStarted(callKey: OutputKey, maybeCallLogs: Option[CallLogs]) extends CallMessage

  sealed trait TerminalCallMessage extends CallMessage

  case class CallAborted(callKey: OutputKey) extends TerminalCallMessage

  case class CallCompleted(callKey: OutputKey, callOutputs: CallOutputs, executionEvents: Seq[ExecutionEventEntry], returnCode: Int, hash: Option[ExecutionHash], resultsClonedFrom: Option[BackendCallJobDescriptor]) extends TerminalCallMessage

  case class ScatterCompleted(callKey: ScatterKey) extends TerminalCallMessage

  case class CallFailedRetryable(callKey: OutputKey, executionEvents: Seq[ExecutionEventEntry], returnCode: Option[Int], failure: Throwable) extends TerminalCallMessage

  case class CallFailedNonRetryable(callKey: OutputKey, executionEvents: Seq[ExecutionEventEntry], returnCode: Option[Int], failure: String) extends TerminalCallMessage

  case class CallFailedToInitialize(callKey: ExecutionStoreKey, reason: String) extends TerminalCallMessage

  case object Terminate extends WorkflowActorMessage

  final case class CachesCreated(startMode: WorkflowStartMessage) extends WorkflowActorMessage

  final case class AsyncFailure(t: Throwable) extends WorkflowActorMessage

  final case class PerformTransition(toState: WorkflowState) extends WorkflowActorMessage

  sealed trait PersistenceMessage extends WorkflowActorMessage {
    def callKey: ExecutionStoreKey

    def executionStatus: ExecutionStatus
  }

  final private[workflow] case class PersistStatus(callKey: OutputKey, status: ExecutionStatus, callOutputs: Option[CallOutputs], returnCode: Option[Int],
                                         hash: Option[ExecutionHash], resultsClonedFrom: Option[BackendCallJobDescriptor], sender: ActorRef, message: TerminalCallMessage) extends WorkflowActorMessage

  final case class PersistenceSucceeded(callKey: ExecutionStoreKey,
                                        executionStatus: ExecutionStatus,
                                        outputs: Option[CallOutputs] = None) extends PersistenceMessage

  final case class PersistenceFailed(callKey: ExecutionStoreKey, executionStatus: ExecutionStatus) extends PersistenceMessage

  /** Used for exploded scatters which create many shards in one shot. */
  final case class PersistencesCompleted(callKeys: Traversable[ExecutionStoreKey],
                                         executionStatus: ExecutionStatus) extends WorkflowActorMessage

  /**
    * This message is sent from the context of the onTransition handler to trigger a `startRunnableCalls` invocation.
    * `startRunnableCalls` determines the calls which are pending persistence, information which is added to the
    * state data and passed forward to subsequent message processing.
    */
  case object StartRunnableCalls extends WorkflowActorMessage

  /**
    * Trigger a Running workflow to check whether or not it is complete (either successfully or not).
    */
  case object CheckForWorkflowComplete extends WorkflowActorMessage

  /**
    * Message sent to self to actually start a call.  This is necessary to synchronize access to the symbol cache
    * during call input retrieval.  The symbol cache is mutable state which should only be accessed in message
    * processing threads.
    * For initial start, assumes an execution is already persisted in Starting.  A restarted/resumed call should
    * be in Running, but the message handler will freshly persist the execution to Starting to side effect writing
    * a new call start time. */
  sealed trait CallStartMessage extends WorkflowActorMessage {
    def callKey: CallKey

    def startMode: CallActor.StartMode
  }

  /** Represents starting a call for the first time, as opposed to a restart. */
  final case class InitialStartCall(override val callKey: CallKey,
                                    override val startMode: CallActor.StartMode) extends CallStartMessage

  /** This signifies using an existing previously run call to fulfill the results of the callKey. */
  final case class UseCachedCall(override val callKey: BackendCallKey,
                                 override val startMode: CallActor.UseCachedCall) extends CallStartMessage

  /** Represents restarting a call for backends which support restart. */
  final case class RestartCall(override val callKey: CallKey, override val startMode: CallActor.StartMode) extends CallStartMessage

  sealed trait WorkflowStartMessage extends WorkflowActorMessage {
    def replyTo: Option[ActorRef]
  }

  case class StartNewWorkflow(replyTo: Option[ActorRef] = None) extends WorkflowStartMessage

  case object RestartWorkflow extends WorkflowStartMessage {
    override def replyTo: Option[ActorRef] = None
  }
}