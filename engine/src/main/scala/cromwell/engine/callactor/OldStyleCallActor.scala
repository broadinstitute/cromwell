package cromwell.engine.callactor

import akka.actor.FSM.NullFunction
import akka.actor._
import akka.event.Logging
import com.google.api.client.util.ExponentialBackOff
import cromwell.engine._
import cromwell.engine.backend._
import cromwell.engine.callactor.OldStyleCallActor.{CallActorData, CallActorState}
import cromwell.engine.callexecution.OldStyleCallExecutionActor
import cromwell.engine.callexecution.OldStyleCallExecutionActor.CallExecutionActorMessage
import cromwell.engine.workflow.{CallKey, OldStyleWorkflowActor}
import cromwell.logging.WorkflowLogger
import wdl4s.Scope
import wdl4s.values.WdlValue

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.language.postfixOps
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object OldStyleCallActor {

  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  sealed trait CallActorMessage
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  sealed trait StartMode extends CallActorMessage {
    def executionMessage: CallExecutionActorMessage
    def maybeCallLogs: Option[CallLogs]
  }
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case object StartFinalCall extends StartMode {
    override val executionMessage = OldStyleCallExecutionActor.Execute
    override val maybeCallLogs = None
  }
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  final case class StartBackendCall(maybeCallLogs: Option[CallLogs]) extends StartMode {
    override val executionMessage = OldStyleCallExecutionActor.Execute
  }
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  final case class Resume(executionInfos: Map[String, Option[String]]) extends StartMode {
    override val executionMessage = OldStyleCallExecutionActor.Resume(executionInfos)
    override val maybeCallLogs = None
  }
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  final case class UseCachedCall(cachedBackendCall: OldStyleBackendCallJobDescriptor, backendCall: OldStyleBackendCallJobDescriptor,
                                 callLogs: CallLogs) extends StartMode {
    override val executionMessage = OldStyleCallExecutionActor.UseCachedCall(cachedBackendCall)
    override val maybeCallLogs = Option(callLogs)
  }
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  final case class RegisterCallAbortFunction(abortFunction: AbortFunction) extends CallActorMessage
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case object AbortCall extends CallActorMessage
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  final case class ExecutionFinished(call: Scope, executionResult: OldStyleExecutionResult) extends CallActorMessage

  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  sealed trait CallActorState
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case object CallNotStarted extends CallActorState
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case object CallRunningAbortUnavailable extends CallActorState
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case object CallRunningAbortAvailable extends CallActorState
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case object CallRunningAbortRequested extends CallActorState
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case object CallAborting extends CallActorState
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case object CallDone extends CallActorState

  /** The `WorkflowActor` will drop `TerminalCallMessage`s that it is unable to immediately process.  This message
    * represents the acknowledgment from the `WorkflowActor` that the `callMessage` was processed and the status
    * change has been successfully committed to the database. */
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  final case class Ack(callMessage: OldStyleWorkflowActor.TerminalCallMessage) extends CallActorMessage

  /** Message to self to retry sending the `callMessage` to the `WorkflowActor` in the absence of an `Ack`. */
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  final case class Retry(callMessage: OldStyleWorkflowActor.CallMessage) extends CallActorMessage

  /** FSM data class with an optional abort function and exponential backoff. */
  @deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
  case class CallActorData(abortFunction: Option[AbortFunction] = None,
                           backoff: Option[ExponentialBackOff] = None,
                           timer: Option[Cancellable] = None) {

    /**
     * Schedule a single retry of sending the specified message to the specified actor.  Return a copy of this
     * `CallActorData` holding a reference to the timer to allow for its cancellation in the event an `Ack` is
     * received prior to its expiration.
     */
    def copyWithRetry(actor: Actor, callMessage: OldStyleWorkflowActor.CallMessage)(implicit ec: ExecutionContext): CallActorData = {
      val timer = actor.context.system.scheduler.scheduleOnce(backoff.get.nextBackOffMillis().millis) {
        actor.self ! Retry(callMessage)
      }
      this.copy(timer = Option(timer))
    }

    /** Attempt to cancel the last timer if there is one registered. */
    def cancelTimer() = timer foreach { _.cancel() }
  }

  def props(backendCallDescriptor: OldStyleBackendCallJobDescriptor): Props = Props(new OldStyleBackendCallActor(backendCallDescriptor))
  def props(finalCallDescriptor: FinalCallJobDescriptor) = Props(new OldStyleFinalCallActor(finalCallDescriptor))
}

/** Actor to manage the execution of a single call. */
// FIXME It feels like I shouldn't need to restate the type bounds of JobDescriptor's CallKey type variable.
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
trait OldStyleCallActor[D <: OldStyleJobDescriptor[_ <: CallKey]] extends LoggingFSM[CallActorState, CallActorData] with CromwellActor {
  def jobDescriptor: D
  def key = jobDescriptor.key
  def workflowDescriptor = jobDescriptor.workflowDescriptor
  protected def callExecutionActor: ActorRef

  import OldStyleCallActor._

  type CallOutputs = Map[String, WdlValue]

  startWith(CallNotStarted, CallActorData())

  implicit val ec = context.system.dispatcher

  val call = key.scope
  val akkaLogger = Logging(context.system, classOf[OldStyleCallActor[D]])
  val logger = WorkflowLogger(
    "CallActor",
    workflowDescriptor,
    akkaLogger = Option(akkaLogger),
    callTag = Option(key.tag)
  )

  // Called on every state transition.
  onTransition {
    case _ -> CallDone => shutDown()
    case fromState -> toState =>
      // Only log this at debug - these states are never seen or used outside of the CallActor itself.
      logger.debug(s"transitioning from $fromState to $toState.")
  }

  when(CallNotStarted) {
    case Event(startMode: StartMode, _) =>
      // There's no special Retry/Ack handling required for CallStarted message, the WorkflowActor can always
      // handle those immediately.
      context.parent ! OldStyleWorkflowActor.CallStarted(key, startMode.maybeCallLogs)
      callExecutionActor ! startMode.executionMessage
      goto(CallRunningAbortUnavailable)
    case Event(AbortCall, _) => handleFinished(call, OldStyleAbortedExecution)
  }

  when(CallRunningAbortUnavailable) {
    case Event(RegisterCallAbortFunction(newAbortFunction), data) =>
      val updatedData = data.copy(abortFunction = Option(newAbortFunction))
      goto(CallRunningAbortAvailable) using updatedData
    case Event(AbortCall, _) => goto(CallRunningAbortRequested)
  }

  when(CallRunningAbortAvailable) {
    case Event(RegisterCallAbortFunction(newAbortFunction), data) =>
      logger.warn("An existing abort function was overwritten with a new one. This might indicate unexpected state transitions in the backend.")
      val updatedData = data.copy(abortFunction = Option(newAbortFunction))
      goto(CallRunningAbortAvailable) using updatedData
    case Event(AbortCall, data) => tryAbort(data)
  }

  when(CallRunningAbortRequested) {
    case Event(RegisterCallAbortFunction(abortFunction), data) =>
      val updatedData = data.copy(abortFunction = Option(abortFunction))
      tryAbort(updatedData)
  }

  // When in CallAborting, the only message being listened for is the CallComplete message (which is already
  // handled by whenUnhandled.)
  when(CallAborting) { NullFunction }

  when(CallDone) {
    case Event(e, _) =>
      logger.error(s"received unexpected event $e while in state $stateName")
      stay()
  }

  whenUnhandled {
    case Event(ExecutionFinished(finishedCall, executionResult), _) => handleFinished(finishedCall, executionResult)
    case Event(Retry(message), data) =>
      logger.debug(s"CallActor retrying ${message.callKey.tag}")
      // Continue retrying this message until it is Acked.
      context.parent ! message
      val updatedData = data.copyWithRetry(this, message)
      stay() using updatedData
    case Event(Ack(message: OldStyleWorkflowActor.TerminalCallMessage), data) =>
      logger.debug(s"CallActor received ack for ${message.callKey.tag}")
      data.cancelTimer()
      val updatedData = data.copy(backoff = None, timer = None)
      goto(CallDone) using updatedData
    case Event(e, _) =>
      logger.error(s"received unhandled event $e while in state $stateName")
      stay()
  }

  private def tryAbort(data: CallActorData): OldStyleCallActor.this.State = {
    data.abortFunction match {
      case Some(af) =>
        logger.info("Abort function called.")
        af.function()
        goto(CallAborting) using data
      case None =>
        logger.warn("Call abort failed because the provided abort function was null.")
        goto(CallRunningAbortRequested) using data
    }
  }

  private def handleFinished(call: Scope, executionResult: OldStyleExecutionResult): State = {

    def createBackoff: Option[ExponentialBackOff] = Option(
      new ExponentialBackOff.Builder()
        .setInitialIntervalMillis(100)
        .setMaxIntervalMillis(3000)
        .setMaxElapsedTimeMillis(Integer.MAX_VALUE)
        .setMultiplier(1.1)
        .build()
    )

    val message = executionResult match {
      case OldStyleSuccessfulBackendCallExecution(outputs, executionEvents, returnCode, hash, resultsClonedFrom) =>
        OldStyleWorkflowActor.CallCompleted(key, outputs, executionEvents, returnCode, if (workflowDescriptor.writeToCache) Option(hash) else None, resultsClonedFrom)
      case OldStyleSuccessfulFinalCallExecution => OldStyleWorkflowActor.CallCompleted(key, Map.empty, Seq.empty, 0, None, None)
      case OldStyleAbortedExecution => OldStyleWorkflowActor.CallAborted(key)
      case OldStyleRetryableFailedExecution(e, returnCode, events) =>
        logger.error("Failing call with retryable Failure: " + e.getMessage, e)
        OldStyleWorkflowActor.CallFailedRetryable(key, events, returnCode, e)
      case OldStyleNonRetryableFailedExecution(e, returnCode, events) =>
        logger.error("Failing call: " + e.getMessage, e)
        OldStyleWorkflowActor.CallFailedNonRetryable(key, events, returnCode, e.getMessage)
    }

    context.parent ! message
    // The WorkflowActor will drop TerminalCallMessages that it is unable to immediately process.
    // Retry sending this message to the WorkflowActor until it is Acked.
    val updatedData = stateData.copy(backoff = createBackoff).copyWithRetry(this, message)
    stay using updatedData
  }

  private def shutDown(): Unit = {
    logger.debug(s"done, shutting down.")
    context.stop(self)
  }
}
