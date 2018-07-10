package wes2cromwell

import akka.actor.{Actor, ActorLogging, Props}
import akka.http.scaladsl.model._
import akka.stream.scaladsl.Source
import akka.util.ByteString
import cromiam.auth.Collection.LabelsKey
import cromiam.webservice.SubmissionSupport.{WorkflowInputsKey, WorkflowOptionsKey, WorkflowTypeKey, WorkflowTypeVersionKey}
import spray.json.JsObject

import scala.concurrent.ExecutionContext.Implicits.global

// FIXME: stuff in wrong file

final case class WesSubmission(
   workflowParams: String,
   workflowType: String,
   workflowTypeVersion: String,
   tags: Option[String],
   workflowEngineParameters: Option[String],
   workflowUrl: String,
   workflowAttachment: Iterable[String]
) {
  val entity: MessageEntity = {
    val sourcePart = ???
    val typePart = Multipart.FormData.BodyPart(WorkflowTypeKey, HttpEntity(MediaTypes.`application/json`, workflowType))
    val typeVersionPart = Multipart.FormData.BodyPart(WorkflowTypeVersionKey, HttpEntity(MediaTypes.`application/json`, workflowTypeVersion))
    val inputsPart = Multipart.FormData.BodyPart(WorkflowInputsKey, HttpEntity(MediaTypes.`application/json`, workflowParams))
    val optionsPart = workflowEngineParameters map { o => Multipart.FormData.BodyPart(WorkflowOptionsKey, HttpEntity(MediaTypes.`application/json`, o)) }
    val labelsPart = tags map { t => Multipart.FormData.BodyPart(LabelsKey, HttpEntity(MediaTypes.`application/json`, t)) }

    val parts = List(Option(sourcePart), Option(typePart), Option(typeVersionPart), Option(inputsPart), optionsPart, labelsPart).flatten

    Multipart.FormData(parts: _*).toEntity()
  }

  /**
    * Needs to handle a few things:
    *
    * If just workflowUrl is specified, pass it along to Cromwell's workflow_url field. This relies on #3849 but not worth the hassle of working around if that doesn't exist yet
    *
    * If workflowAttachment is populated, then workflowUrl becomes a relative URL into this Iterable. For now,
    * this means taking *that* member and using it as the sourcePart and then zip encoding the rest to go into the zip bundle for Cromwell
    */
  private def getSourceParts(): Unit = {

  }
}



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
