package cromwell.engine

import akka.actor.{LoggingFSM, FSM, Props}
import akka.event.Logging
import akka.pattern.{ask, pipe}
import akka.util.Timeout
import cromwell.binding.values.WdlValue
import cromwell.binding._
import cromwell.engine.StoreActor._
import cromwell.engine.backend.Backend
import WorkflowActor._
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration._
import scala.language.postfixOps

object WorkflowActor {
  sealed trait WorkflowActorMessage
  case object Complete extends WorkflowActorMessage
  case object Start extends WorkflowActorMessage
  case class Failed(failure: String) extends WorkflowActorMessage
  case object GetOutputs extends WorkflowActorMessage
  case object GetFailureMessage extends WorkflowActorMessage
  

  def props(id: WorkflowId, namespace: WdlNamespace, coercedInputs: WorkflowCoercedInputs, backend: Backend): Props = {
    Props(new WorkflowActor(id, namespace, coercedInputs, backend))
  }

  implicit val ActorTimeout = Timeout(5 seconds)

  sealed trait WorkflowFailure
  case object NoFailureMessage extends WorkflowFailure
  case class FailureMessage(msg: String) extends WorkflowFailure with WorkflowActorMessage
}

case class WorkflowActor(id: WorkflowId,
                         namespace: WdlNamespace,
                         actualInputs: WorkflowCoercedInputs,
                         backend: Backend) extends LoggingFSM[WorkflowState, WorkflowFailure] {
  implicit val timeout = Timeout(5 seconds)

  private val storeActor = context.actorOf(StoreActor.props(namespace, actualInputs))

  startWith(WorkflowSubmitted, NoFailureMessage)

  when(WorkflowSubmitted) {
    case Event(WorkflowActor.Start, NoFailureMessage) => goto(WorkflowRunning)
  }

  when(WorkflowRunning) {
    case Event(CallActor.Started(call), NoFailureMessage) =>
      storeActor ! UpdateStatus(call, ExecutionStatus.Running)
      stay()
    case Event(CallActor.Completed(completedCall, callOutputs), NoFailureMessage) =>
      storeActor ! CallCompleted(completedCall, callOutputs)
      stay()
    case Event(StoreActor.RunnableCalls(runnableCalls), NoFailureMessage) =>
      if (runnableCalls.nonEmpty) {
        log.info("Starting calls: " + runnableCalls.map {_.name}.toSeq.sorted.mkString(", "))
      }
      runnableCalls foreach startCallActor
      stay()
    case Event(CallActor.Failed(call, failure), NoFailureMessage) =>
      storeActor ! UpdateStatus(call, ExecutionStatus.Failed)
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
    case WorkflowSubmitted -> WorkflowRunning => storeActor ! FindRunnableCalls
  }

  /** Create a per-call `CallActor` for the specified `Call` and send it a `Start` message to
    * begin execution. */
  private def startCallActor(call: Call): Unit = {
    val callActorProps = CallActor.props(call, backend, storeActor, "CallActor-" + call.name)
    context.actorOf(callActorProps) ! CallActor.Start
  }
}