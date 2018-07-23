package wes2cromwell

import java.net.URL

import akka.actor.{Actor, ActorLogging, Props}

import scala.concurrent.ExecutionContext.Implicits.global

// FIXME: remove all of this

object WorkflowActor {
  final case object GetWorkflows
  final case class GetWorkflow(workflowId: String)

  def props: Props = Props[WorkflowActor]
}

class WorkflowActor extends Actor with ActorLogging {
  import WorkflowActor._
  lazy val transmogriphy = new Wes2CromwellInterface(new URL("http://some.bullshit.com"))(context.system, global)

  def receive: Receive = {
    case GetWorkflows =>
      transmogriphy.getWorkflows(sender())
  }
}
