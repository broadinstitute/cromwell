package cromwell.engine

import akka.actor.Actor
import akka.actor.SupervisorStrategy.Stop

/**
 * This model assumes workflow actors are transient, created once per workflow instance.
 */
object WorkflowActor {

  // This is hand-waving over the actual engine-level representation of a workflow.
  type WorkflowModel = String

  case class Construct(model: WorkflowModel)
  case object Constructed

  case object Start
  case object Started

  // Note there is no Stop message defined here, Stop is from Akka.
  case object Stopped
}

// Represents the root of a single workflow instance, not a manager of multiple
// workflows.
class WorkflowActor extends Actor {
  import WorkflowActor._

  override def receive: Receive = {

    case Construct(model) =>
      sender ! Constructed

    case Start =>
      sender ! Started

    case Stop =>
      // TODO is this the right way to stop this actor?
      // TODO http://stackoverflow.com/questions/13847963/akka-kill-vs-stop-vs-poison-pill
      context.stop(self)
      sender ! Stopped
  }
}
