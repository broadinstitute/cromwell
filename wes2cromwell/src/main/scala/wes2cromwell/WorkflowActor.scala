package wes2cromwell

import akka.actor.{ Actor, ActorLogging, Props }
import spray.json.JsObject

import scala.concurrent.ExecutionContext.Implicits.global

final case class WorkflowRequest(
  workflow_descriptor: Option[String],
  workflow_params: Option[JsObject],
  workflow_type: String,
  workflow_type_version: String,
  tags: Option[JsObject],
  workflow_engine_parameters: Option[JsObject],
  workflow_url: Option[String]
)

final case class WorkflowDescription(
  workflow_id: String,
  state: WorkflowState
)

final case class WorkflowTypeVersion(workflow_type_version: Seq[String])

final case class ErrorResponse(
  msg: String,
  status_code: Int
)

object WorkflowActor {
  final case object GetWorkflows
  final case class PostWorkflow(workflowRequest: WorkflowRequest)
  final case class GetWorkflow(workflowId: String)
  final case class DeleteWorkflow(workflowId: String)
  final case class GetWorkflowStatus(workflowId: String)

  def props: Props = Props[WorkflowActor]
}

class WorkflowActor extends Actor with ActorLogging {
  import WorkflowActor._
  lazy val transmogriphy = new Transmogriphy()(context.system, global)

  def receive: Receive = {
    case GetWorkflows =>
      transmogriphy.getWorkflows(sender())
    case PostWorkflow(workflowRequest) =>
      transmogriphy.postWorkflow(sender(), workflowRequest)
    case GetWorkflow(workflowId) =>
      transmogriphy.getWorkflow(sender(), workflowId)
    case DeleteWorkflow(workflowId) =>
      transmogriphy.deleteWorkflow(sender(), workflowId)
    case GetWorkflowStatus(workflowId) =>
      transmogriphy.getWorkflowStatus(sender(), workflowId)
  }
}
