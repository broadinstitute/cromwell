package cromwell.engine

import akka.actor.FSM.NullFunction
import akka.actor.{LoggingFSM, Props}
import cromwell.binding._
import cromwell.binding.values.WdlValue
import cromwell.engine.CallActor.CallActorState
import cromwell.engine.backend.{Backend, TaskAbortedException}
import cromwell.engine.workflow.WorkflowActor.CallFailed
import cromwell.engine.workflow.{CallKey, WorkflowActor}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object CallActor {

  sealed trait CallActorMessage
  case object Start extends CallActorMessage
  final case class RegisterCallAbortFunction(abortFunction: AbortFunction) extends CallActorMessage
  case object AbortCall extends CallActorMessage
  final case class ExecutionFinished(call: Call, outputs: Try[CallOutputs]) extends CallActorMessage

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
  val callReference = CallReference(workflowDescriptor.name, workflowDescriptor.id, call.fullyQualifiedName)
  val tag = s"CallActor [$callReference]"

  // Called on every state transition.
  onTransition {
    case fromState -> toState =>
      // Only log this at debug - these states are never seen or used outside of the CallActor itself.
      log.debug(s"$tag transitioning from $fromState to $toState.")
  }

  when(CallNotStarted) {
    case Event(Start, _) =>
      sender() ! WorkflowActor.CallStarted(key)
      val backendCall = backend.bindCall(workflowDescriptor, key, locallyQualifiedInputs, AbortRegistrationFunction(registerAbortFunction(callReference)))
      val executionActorName = s"CallExecutionActor-${workflowDescriptor.id}-${call.name}"
      context.actorOf(CallExecutionActor.props(backendCall), executionActorName) ! CallExecutionActor.Execute
      goto(CallRunningAbortUnavailable)
    case Event(AbortCall, _) => handleFinished(call, Failure(new TaskAbortedException()))
  }

  when(CallRunningAbortUnavailable) {
    case Event(RegisterCallAbortFunction(abortFunction), _) => goto(CallRunningAbortAvailable) using Option(abortFunction)
    case Event(AbortCall, _) => goto(CallRunningAbortRequested)
  }

  when(CallRunningAbortAvailable) {
    case Event(AbortCall, abortFunction) => tryAbort(abortFunction)
  }

  when(CallRunningAbortRequested) {
    case Event(RegisterCallAbortFunction(abortFunction: AbortFunction), _) => tryAbort(Option(abortFunction))
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
    case Event(ExecutionFinished(call: Call, outputs: Try[CallOutputs]), _) => handleFinished(call, outputs)
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

  private def handleFinished(call: Call, outputTry: Try[CallOutputs]): CallActor.this.State = {
    outputTry match {
      case Success(outputs) => context.parent ! WorkflowActor.CallCompleted(key, outputs)
      case Failure(e: TaskAbortedException) =>
            context.parent ! WorkflowActor.AbortComplete(key)
      case Failure(e: Throwable) =>
            context.parent ! WorkflowActor.CallFailed(key, e.getMessage)
    }

    goto(CallDone)
  }

  private def registerAbortFunction(callReference: CallReference)(abortFunction: AbortFunction): Unit = {
    self ! CallActor.RegisterCallAbortFunction(abortFunction)
  }
}
