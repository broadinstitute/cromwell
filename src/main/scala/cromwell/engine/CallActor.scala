package cromwell.engine

import akka.actor.FSM.NullFunction
import akka.actor.{LoggingFSM, Props}
import cromwell.binding._
import cromwell.binding.values.WdlValue
import cromwell.engine.CallActor.CallActorState
import cromwell.engine.backend.{TaskAbortedException, Backend}
import cromwell.engine.workflow.WorkflowActor
import cromwell.engine.workflow.WorkflowActor.CallFailed

import scala.language.postfixOps
import scala.util.{Try, Failure, Success}

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

  def props(call: Call, locallyQualifiedInputs: CallInputs, backend: Backend, workflowDescriptor: WorkflowDescriptor): Props =
    Props(new CallActor(call, locallyQualifiedInputs, backend, workflowDescriptor))
}

/** Actor to manage the execution of a single call. */
class CallActor(call: Call, locallyQualifiedInputs: CallInputs, backend: Backend, workflowDescriptor: WorkflowDescriptor)
  extends LoggingFSM[CallActorState, Option[AbortFunction]] with CromwellActor {

  import CallActor._

  type CallOutputs = Map[String, WdlValue]

  startWith(CallNotStarted, None)

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
      val backendInputs = backend.adjustInputPaths(call, locallyQualifiedInputs)
      call.instantiateCommandLine(backendInputs, backend.setupCallEnvironment(call, workflowDescriptor).engineFunctions) match {
        case Success(commandLine) =>
          sender() ! WorkflowActor.CallStarted(call)
          context.actorOf(CallExecutionActor.props(callReference)) ! CallExecutionActor.Execute(workflowDescriptor.id, backend, commandLine, workflowDescriptor, call, backendInputs, inputName => locallyQualifiedInputs.get(inputName).get)
          goto(CallRunningAbortUnavailable)
        case Failure(e) =>
          val message = s"Call '${call.fullyQualifiedName}' failed to launch command: " + e.getMessage
          log.error(e, s"$tag: $call failed: $message")
          context.parent ! CallFailed(call, message)
          goto(CallDone)
      }
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
      case Success(outputs) => context.parent ! WorkflowActor.CallCompleted(call, outputs)
      case Failure(e: TaskAbortedException) =>
            context.parent ! WorkflowActor.AbortComplete(call)
      case Failure(e: Throwable) =>
            context.parent ! WorkflowActor.CallFailed(call, e.getMessage)
    }

    goto(CallDone)
  }
}
