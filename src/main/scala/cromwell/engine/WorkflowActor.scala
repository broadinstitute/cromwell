package cromwell.engine

import akka.actor.SupervisorStrategy.Stop
import akka.actor.{Actor, ActorRef, Props}
import akka.event.LoggingReceive
import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, FullyQualifiedName, WdlBinding}

/**
 * This model assumes workflow actors are transient, created once per workflow instance.
 */
object WorkflowActor {

  sealed trait WorkflowActorMessage

  case object Start extends WorkflowActorMessage

  case object Started extends WorkflowActorMessage

  // Note there is no Stop message defined here, Stop is from Akka.
  case object Stopped extends WorkflowActorMessage

  case class Done(outputs: Map[String, WdlValue]) extends WorkflowActorMessage

  case class InvalidOperation(message: String) extends WorkflowActorMessage

  case class Failed(throwable: Throwable) extends WorkflowActorMessage

  def buildWorkflowActorProps(binding: WdlBinding, actualInputs: Map[FullyQualifiedName, WdlValue]): Props = {
    // Check inputs for presence and type compatibility
    val diagnostics = binding.workflow.inputs.collect {
      case requiredInput if !actualInputs.contains(requiredInput._1) =>
        requiredInput._1 -> "Required workflow input not specified"

      case requiredInput if actualInputs.get(requiredInput._1).get.wdlType != requiredInput._2 =>
        val expected = actualInputs.get(requiredInput._1).get.wdlType
        // FIXME formatting
        requiredInput._1 -> s"Incompatible workflow input types, expected $expected, got ${requiredInput._2}"
    }

    if (diagnostics.nonEmpty) throw new UnsatisfiedInputsException(diagnostics)

    Props(new WorkflowActor(binding, actualInputs))
  }
}

/** Represents the root of a single workflow instance, not a manager of multiple
  * workflows.
  */
class WorkflowActor private(binding: WdlBinding, actualInputs: Map[FullyQualifiedName, WdlValue]) extends Actor {

  import WorkflowActor._

  // That which started the workflow, a role which has not yet been defined more specifically.
  private var _primeMover: Option[ActorRef] = None

  private def isStarted: Boolean = _primeMover.isDefined

  private def primeMover = _primeMover.get

  override def receive: Receive = {
    LoggingReceive {

      case Start =>
        _primeMover = Option(sender())
        val executionStatusStore = new ExecutionStatusStore(binding)
        val symbolStore = new SymbolStore(binding, actualInputs)
        executionStatusStore.runnableCalls.foreach { call =>
          startCallActor(symbolStore, call)
        }
        primeMover ! Started

      case CallActor.Started =>
      // TODO update execution status store

      case CallActor.Done(taskOutputs) =>
        // TODO Update execution status and symbol stores, evaluate workflow outputs
        // TODO There is only one runnable task at the moment, so immediately
        // TODO report as done to the creator.
        primeMover ! Done(Map.empty[String, WdlValue])
      // TODO This should reevaluate runnable tasks, if there are no more
      // TODO initiate shutdown, including waiting for child Call actors to
      // TODO report as Stopped (state "SHUTTING_DOWN"?).  Shutting this down
      // TODO before waiting for child actors to report as stopped
      // TODO would produce a slew of dead letter warnings as a shutdown
      // TODO of self would terminate the entire actor hierarchy.
      // TODO If there are more runnable tasks they should be started and
      // TODO the stores updated appropriately.
      // TODO The code needed here is not very different from what Start
      // TODO should be doing.

      case CallActor.Failed(t) =>
        // TODO update execution status store
        primeMover ! Failed(t)

      case CallActor.Stopped =>
      // TODO update execution status store

      case Stop if !isStarted =>
        sender ! InvalidOperation("Stop issued against WorkflowActor that appears never to have been started")

      case Stop =>
        // TODO is this the right way to stop this actor?
        // TODO http://stackoverflow.com/questions/13847963/akka-kill-vs-stop-vs-poison-pill
        context.stop(self)
        primeMover ! Stopped

      case unknown@_ =>
        primeMover ! InvalidOperation(s"Unknown message '$unknown'")
    }
  }

  private def startCallActor(symbolStore: SymbolStore, call: Call): Unit = {
    val callActor = context.actorOf(CallActor.props, "CallActor-" + call.name)
    callActor ! CallActor.Start(call, symbolStore)
  }
}
