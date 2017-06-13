package cromwell.webservice

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.labels.Labels
import cromwell.core.{WorkflowId, WorkflowMetadataKeys}
import cromwell.services.metadata.MetadataService._
import cromwell.webservice.LabelManagerActor._
import cromwell.webservice.PerRequest.RequestComplete
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._
import spray.json.{DefaultJsonProtocol, JsObject, JsString}


object LabelManagerActor {

  def props(serviceRegistryActor: ActorRef) = Props(new LabelManagerActor(serviceRegistryActor))

  final case class LabelData(workflowId: WorkflowId, labels: Labels)

  sealed trait LabelMessage {
    def data: LabelData
  }

  sealed trait LabelAction extends LabelMessage
  final case class AddLabel(data: LabelData) extends LabelAction

  sealed trait LabelResponse extends LabelMessage

  def processLabelsResponse(workflowId: String, labels: Map[String, String]): JsObject = {
    JsObject(Map(
      WorkflowMetadataKeys.Id -> JsString(workflowId.toString),
      WorkflowMetadataKeys.Labels -> JsObject(labels map { case (k, v) => k -> JsString(v) })
    ))
  }
}

class LabelManagerActor(serviceRegistryActor: ActorRef) extends Actor with ActorLogging with DefaultJsonProtocol {

  implicit val ec = context.dispatcher

  import WorkflowJsonSupport._

  def receive = {
    case add: AddLabelsToWorkflowMetadata =>
      serviceRegistryActor ! add
    case LabelUpdateSuccess(id, labels) =>
      val response = processLabelsResponse(id, labels)
      context.parent ! RequestComplete((StatusCodes.OK, response))
    case failure: LabelUpdateFailure =>
      val response = APIResponse.fail(new RuntimeException(s"Can't find metadata service to update labels for ${failure.id} due to ${failure.reason.getMessage}"))
      context.parent ! RequestComplete((StatusCodes.InternalServerError, response))
    case random => log.error("Unexpected message to LabelManagerActor: {}", random)
  }

}
