package cromwell.engine

import akka.actor.{FSM, LoggingFSM, Props}
import akka.pattern.{ask, pipe}
import cromwell.binding._
import cromwell.engine.WorkflowActor._
import cromwell.engine.backend.Backend

import scala.concurrent.ExecutionContext.Implicits.global
import scala.language.postfixOps

object WorkflowActor {
  sealed trait WorkflowActorMessage
  case object Start extends WorkflowActorMessage
  case object Complete extends WorkflowActorMessage
  case object GetFailureMessage extends WorkflowActorMessage
  case object GetOutputs extends WorkflowActorMessage
  case class CallStarted(call: Call) extends WorkflowActorMessage
  case class CallCompleted(call: Call, outputs: WorkflowOutputs) extends WorkflowActorMessage
  case class CallFailed(call: Call, failure: String) extends WorkflowActorMessage
  case class RunnableCalls(calls: Iterable[Call]) extends WorkflowActorMessage

  def props(id: WorkflowId, namespace: WdlNamespace, coercedInputs: WorkflowCoercedInputs, backend: Backend): Props = {
    Props(new WorkflowActor(id, namespace, coercedInputs, backend))
  }

  sealed trait WorkflowFailure
  case object NoFailureMessage extends WorkflowFailure
  case class FailureMessage(msg: String) extends WorkflowFailure with WorkflowActorMessage
}

case class WorkflowActor(id: WorkflowId,
                         namespace: WdlNamespace,
                         actualInputs: WorkflowCoercedInputs,
                         backend: Backend) extends LoggingFSM[WorkflowState, WorkflowFailure] with CromwellActor {
  private val storeActor = context.actorOf(StoreActor.props(namespace, actualInputs))

  startWith(WorkflowSubmitted, NoFailureMessage)

  when(WorkflowSubmitted) {
    case Event(WorkflowActor.Start, NoFailureMessage) => goto(WorkflowRunning)
  }

  when(WorkflowRunning) {
    case Event(CallStarted(call), NoFailureMessage) =>
      storeActor ! StoreActor.UpdateStatus(call, ExecutionStatus.Running)
      stay()
    case Event(CallCompleted(completedCall, callOutputs), NoFailureMessage) =>
      storeActor ! StoreActor.CallCompleted(completedCall, callOutputs)
      stay()
    case Event(RunnableCalls(runnableCalls), NoFailureMessage) =>
      if (runnableCalls.nonEmpty) {
        log.info("Starting calls: " + runnableCalls.map {_.name}.toSeq.sorted.mkString(", "))
      }
      runnableCalls foreach startCallActor
      stay()
    case Event(CallFailed(call, failure), NoFailureMessage) =>
      storeActor ! StoreActor.UpdateStatus(call, ExecutionStatus.Failed)
      goto(WorkflowFailed) using FailureMessage(failure)
    case Event(WorkflowActor.Complete, NoFailureMessage) => goto(WorkflowSucceeded)
  }

  when(WorkflowFailed) {
    case Event(GetFailureMessage, msg: FailureMessage) =>
      sender() ! msg
      stay()
  }

  // We're supporting GetOutputs for all states, so there's nothing particular to do here
  when(WorkflowSucceeded)(FSM.NullFunction)

  whenUnhandled {
    case Event(WorkflowActor.GetOutputs, _) =>
      storeActor.ask(StoreActor.GetOutputs) pipeTo sender
      stay()
    case Event(e, _) =>
      log.debug(s"Received unhandled event $e while in state $stateName")
      stay()
  }

  onTransition {
    case WorkflowSubmitted -> WorkflowRunning => storeActor ! StoreActor.StartRunnableCalls
  }

  /** Create a per-call `CallActor` for the specified `Call` and send it a `Start` message to
    * begin execution. */
  private def startCallActor(call: Call): Unit = {
    val callActorProps = CallActor.props(call, backend, storeActor, "CallActor-" + call.name)
    context.actorOf(callActorProps) ! CallActor.Start
  }
}