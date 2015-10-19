package cromwell.engine

import akka.actor.FSM.NullFunction
import akka.actor.{LoggingFSM, Props}
import cromwell.binding._
import cromwell.binding.values.WdlValue
import cromwell.engine.CallActor.CallActorState
import cromwell.engine.CallExecutionActor.CallExecutionActorMessage
import cromwell.engine.backend._
import cromwell.engine.workflow.{CallKey, WorkflowActor}

import scala.language.postfixOps

object CallActor {

  sealed trait CallActorMessage
  sealed trait StartMode extends CallActorMessage {
    def executionMessage: CallExecutionActorMessage
  }
  case object Start extends StartMode { override val executionMessage = CallExecutionActor.Execute }
  final case class Resume(jobKey: JobKey) extends StartMode { override val executionMessage = CallExecutionActor.Resume(jobKey) }
  final case class RegisterCallAbortFunction(abortFunction: AbortFunction) extends CallActorMessage
  case object AbortCall extends CallActorMessage
  final case class ExecutionFinished(call: Call, executionResult: ExecutionResult) extends CallActorMessage

  sealed trait CallActorState
  case object CallNotStarted extends CallActorState
  case object CallRunningAbortUnavailable extends CallActorState
  case object CallRunningAbortAvailable extends CallActorState
  case object CallRunningAbortRequested extends CallActorState
  case object CallAborting extends CallActorState
  case object CallDone extends CallActorState

  def props(key: CallKey, locallyQualifiedInputs: CallInputs, backend: Backend, workflowDescriptor: WorkflowDescriptor): Props =
    Props(new CallActor(key, locallyQualifiedInputs, backend, workflowDescriptor))
}

/** Actor to manage the execution of a single call. */
class CallActor(key: CallKey, locallyQualifiedInputs: CallInputs, backend: Backend, workflowDescriptor: WorkflowDescriptor)
  extends LoggingFSM[CallActorState, Option[AbortFunction]] with CromwellActor {

  import CallActor._

  type CallOutputs = Map[String, WdlValue]

  startWith(CallNotStarted, None)

  val call = key.scope
  val tag = s"CallActor [UUID(${workflowDescriptor.shortId}):${key.tag}]"

  // Called on every state transition.
  onTransition {
    case _ -> CallDone =>
      log.debug(s"$tag CallActor is done, shutting down.")
      context.stop(self)
    case fromState -> toState =>
      // Only log this at debug - these states are never seen or used outside of the CallActor itself.
      log.debug(s"$tag transitioning from $fromState to $toState.")
  }

  when(CallNotStarted) {
    case Event(startMode: StartMode, _) =>
      sender() ! WorkflowActor.CallStarted(key)
      val backendCall = backend.bindCall(workflowDescriptor, key, locallyQualifiedInputs, AbortRegistrationFunction(registerAbortFunction))
      val executionActorName = s"CallExecutionActor-${workflowDescriptor.id}-${call.name}"
      context.actorOf(CallExecutionActor.props(backendCall), executionActorName) ! startMode.executionMessage
      goto(CallRunningAbortUnavailable)
    case Event(AbortCall, _) => handleFinished(call, AbortedExecution)
  }

  when(CallRunningAbortUnavailable) {
    case Event(RegisterCallAbortFunction(abortFunction), _) => goto(CallRunningAbortAvailable) using Option(abortFunction)
    case Event(AbortCall, _) => goto(CallRunningAbortRequested)
  }

  when(CallRunningAbortAvailable) {
    case Event(AbortCall, abortFunction) => tryAbort(abortFunction)
  }

  when(CallRunningAbortRequested) {
    case Event(RegisterCallAbortFunction(abortFunction), _) => tryAbort(Option(abortFunction))
  }

  // When in CallAborting, the only message being listened for is the CallComplete message (which is already
  // handled by whenUnhandled.)
  when(CallAborting) { NullFunction }

  when(CallDone) {
    case Event(e, _) =>
      log.warning(s"$tag received unexpected event $e while in state $stateName")
      stay()
  }

  whenUnhandled {
    case Event(ExecutionFinished(finishedCall, executionResult), _) => handleFinished(finishedCall, executionResult)
    case Event(e, _) =>
      log.warning(s"$tag received unhandled event $e while in state $stateName")
      stay()
  }

  private def tryAbort(abortFunction: Option[AbortFunction]): CallActor.this.State = {
    abortFunction match {
      case Some(af) =>
        log.info("Abort function called.")
        af.function()
        goto(CallAborting) using abortFunction
      case None =>
        log.warning("Call abort failed because the provided abort function was null.")
        goto(CallRunningAbortRequested) using abortFunction
    }
  }

  private def handleFinished(call: Call, executionResult: ExecutionResult): CallActor.this.State = {
    executionResult match {
      case SuccessfulExecution(outputs, returnCode) => context.parent ! WorkflowActor.CallCompleted(key, outputs, returnCode)
      case AbortedExecution => context.parent ! WorkflowActor.AbortComplete(key)
      case FailedExecution(e, returnCode) => context.parent ! WorkflowActor.CallFailed(key, returnCode, e.getMessage)
    }

    goto(CallDone)
  }

  private def registerAbortFunction(abortFunction: AbortFunction): Unit = {
    self ! CallActor.RegisterCallAbortFunction(abortFunction)
  }
}
