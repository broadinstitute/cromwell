package cromwell.engine

import akka.actor.{FSM, LoggingFSM, Props}
import akka.event.Logging
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

  def props(descriptor: WorkflowDescriptor, backend: Backend) = Props(new WorkflowActor(descriptor, backend))

  sealed trait WorkflowFailure
  case object NoFailureMessage extends WorkflowFailure
  case class FailureMessage(msg: String) extends WorkflowFailure with WorkflowActorMessage
}

case class WorkflowActor(workflow: WorkflowDescriptor, backend: Backend) extends FSM[WorkflowState, WorkflowFailure] with CromwellActor {

  private val storeActor = context.actorOf(StoreActor.props(workflow, backend.initializeForWorkflow(workflow)))
  val tag: String = s"WorkflowActor [UUID(${workflow.shortId})]"
  override val log = Logging(context.system, classOf[WorkflowActor])

  startWith(WorkflowSubmitted, NoFailureMessage)

  when(WorkflowSubmitted) {
    case Event(Start, NoFailureMessage) => goto(WorkflowRunning)
  }

  when(WorkflowRunning) {
    case Event(CallStarted(call), NoFailureMessage) => updateCallStatusToRunning(call)
    case Event(CallCompleted(completedCall, callOutputs), NoFailureMessage) => updateCallStatusToCompleted(completedCall, callOutputs)
    case Event(RunnableCalls(runnableCalls), NoFailureMessage) => receiveRunnableCalls(runnableCalls)
    case Event(CallFailed(call, failure), NoFailureMessage) => updateCallStatusToFailed(call, failure)
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
    case Event(GetOutputs, _) =>
      storeActor.ask(StoreActor.GetOutputs) pipeTo sender
      stay()
    case Event(e, _) =>
      log.debug(s"Received unhandled event $e while in state $stateName")
      stay()
  }

  onTransition {
    case WorkflowSubmitted -> WorkflowRunning => storeActor ! StoreActor.StartRunnableCalls
  }

  private def updateCallStatusToRunning(call: Call) = {
    log.info(s"$tag: call '${call.name}' running")
    storeActor ! StoreActor.UpdateStatus(call, ExecutionStatus.Running)
    stay()
  }

  private def updateCallStatusToCompleted(call: Call, outputs: WorkflowOutputs) = {
    log.info(s"$tag: call '${call.name}' completed")
    storeActor ! StoreActor.CallCompleted(call, outputs)
    stay()
  }

  private def updateCallStatusToFailed(call: Call, failure: String) = {
    log.info(s"$tag: call '${call.name}' failed ($failure)")
    storeActor ! StoreActor.UpdateStatus(call, ExecutionStatus.Failed)
    goto(WorkflowFailed) using FailureMessage(failure)
  }

  private def unknownMessage(e: Any) = {
    log.warning(s"$tag: Unexpected message: $e")
    stay()
  }

  private def receiveRunnableCalls(calls: Iterable[Call]) = {
    if (calls.nonEmpty) {
      log.info(s"$tag: starting " + calls.map {_.name}.toSeq.sorted.mkString(", "))
    }
    calls foreach { call =>
      log.info(s"$tag: launching CallActor for '${call.name}'")
      val callActorProps = CallActor.props(call, backend, workflow, storeActor)
      context.actorOf(callActorProps) ! CallActor.Start
    }
    stay()
  }
}
