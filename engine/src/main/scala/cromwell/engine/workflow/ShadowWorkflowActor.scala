package cromwell.engine.workflow

import akka.actor.{Props, ActorLogging, Actor}

object ShadowWorkflowActor {
  def props(): Props = Props(new ShadowWorkflowActor)
}

/**
  * Class that will [eventually] orchestrate a single workflow.
  */
class ShadowWorkflowActor extends Actor with ActorLogging{

  private val tag = this.getClass.getSimpleName

  override def receive: Receive = {
    case unhandledMessage: Any => log.info(s"$tag received an unhandled message: $unhandledMessage")
  }
}
