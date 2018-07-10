package wes2cromwell

import akka.actor.{ Actor, ActorLogging, Props }
import spray.json.JsObject

import scala.concurrent.ExecutionContext.Implicits.global

// FIXME: stuff in wrong file

final case class WesSubmission(
   workflowParams: Option[JsObject],
   workflowType: String,
   workflowTypeVersion: String,
   tags: Option[JsObject],
   workflowEngineParameters: Option[JsObject],
   workflowUrl: String,
   workflowAttachment: List[String]
)

// FIXME: Might not be used
final case class WorkflowDescription(
  workflow_id: String,
  state: WorkflowState
)

// FIXME: Might not be used
final case class WorkflowTypeVersion(workflow_type_version: Seq[String])

// FIXME: Might not be used
final case class ErrorResponse(
  msg: String,
  status_code: Int
)

object WorkflowActor {
  final case object GetWorkflows
  final case class PostWorkflow(wesSubmission: WesSubmission)
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
