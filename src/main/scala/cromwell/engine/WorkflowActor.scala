package cromwell.engine

import akka.actor.Actor
import akka.actor.SupervisorStrategy.Stop
import cromwell.binding.{Workflow, WdlBinding, WdlValue}

/**
 * This model assumes workflow actors are transient, created once per workflow instance.
 */
object WorkflowActor {

  case class Construct(binding: WdlBinding, inputs: Map[String, WdlValue])
  case object Constructed

  /**
   * @param badInputs A `Map` of locally qualified input names to a corresponding diagnostic message.
   */
  case class ConstructionFailed(badInputs: Map[String, String])

  case object Start
  case object Started

  // Note there is no Stop message defined here, Stop is from Akka.
  case object Stopped

  case object Done

  case class InvalidOperation(message: String)
}

/** Represents the root of a single workflow instance, not a manager of multiple
  * workflows.
  */
class WorkflowActor extends Actor {

  import WorkflowActor._

  // This var would be unnecessary in a FSM that #becomes constructed.
  private var _constructed = false

  /**
   * Check inputs for presence and type compatibility.
   * @param workflow Workflow
   * @param actualInputs Map of <b>local</b> input names to `WdlValue`s
   * @return Map of locally qualified input names to a corresponding diagnostic message.
   */
  def checkUnsatisfiedInputs(workflow: Workflow, actualInputs: Map[String, WdlValue]): Map[String, String] =
    workflow.inputs.collect {

      case requiredInput if !actualInputs.contains(requiredInput._1) =>
        requiredInput._1 -> "Required workflow input not specified"

      case requiredInput if actualInputs.get(requiredInput._1).get.wdlType != requiredInput._2 =>
        val expected = actualInputs.get(requiredInput._1).get.wdlType
        // FIXME formatting
        requiredInput._1 -> s"Incompatible workflow input types, expected $expected, got ${requiredInput._2}"
    }

  override def receive: Receive = {

    case Construct(binding, actualInputs) if !_constructed =>
      val unsatisfiedInputs = checkUnsatisfiedInputs(binding.workflow, actualInputs)
      if (unsatisfiedInputs.isEmpty) {
        _constructed = true
        sender ! Constructed
      } else {
        sender ! ConstructionFailed(unsatisfiedInputs)
      }

    case msg if !_constructed =>
      sender ! InvalidOperation(s"Received $msg but WorkflowActor not Constructed!")

    case msg@Construct(_, _) if _constructed =>
      sender ! InvalidOperation(s"Received $msg but WorkflowActor already Constructed!")

    case Start =>
      sender ! Started

    case Stop =>
      // TODO is this the right way to stop this actor?
      // TODO http://stackoverflow.com/questions/13847963/akka-kill-vs-stop-vs-poison-pill
      context.stop(self)
      sender ! Stopped
  }
}
