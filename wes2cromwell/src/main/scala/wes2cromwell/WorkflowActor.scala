package wes2cromwell

import java.net.URL

import akka.actor.{Actor, ActorLogging, Props}

import scala.concurrent.ExecutionContext.Implicits.global

final case class WorkflowDescription(
  workflow_id: String,
  state: WesRunState
)

final case class WorkflowTypeVersion(workflow_type_version: Seq[String])

final case class ErrorResponse(
  msg: String,
  status_code: Int
)

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
    case GetWorkflow(workflowId) =>
      transmogriphy.getWorkflow(sender(), workflowId)

  }
}
