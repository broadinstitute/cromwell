package cromwell.webservice

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.labels.Labels
import cromwell.core.{WorkflowId, WorkflowMetadataKeys}
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.services.metadata.MetadataService._
import cromwell.webservice.LabelManagerActor._
import cromwell.webservice.PerRequest.RequestComplete
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._
import spray.json.{DefaultJsonProtocol, JsObject, JsString}

import scala.language.postfixOps


object LabelManagerActor {

  def props(serviceRegistryActor: ActorRef) = Props(new LabelManagerActor(serviceRegistryActor))

  final case class LabelData(workflowId: WorkflowId, labels: Labels)

  sealed trait LabelMessage {
    def data: LabelData
  }

  sealed trait LabelAction extends LabelMessage
  final case class LabelAddition(data: LabelData) extends LabelAction

  sealed trait LabelResponse extends LabelMessage

  def processLabelsResponse(workflowId: WorkflowId, labels: Map[String, String]): JsObject = {
    JsObject(Map(
      WorkflowMetadataKeys.Id -> JsString(workflowId.toString),
      WorkflowMetadataKeys.Labels -> JsObject(labels mapValues JsString.apply)
    ))
  }

  def metadataEventsToLabels(events: Iterable[MetadataEvent]): Map[String, String] = {
    events map { x => x.key.key.split("\\:").last -> x.value.get.toString } toMap
  }

  def labelsToMetadataEvents(labels: Labels, workflowId: WorkflowId): Iterable[MetadataEvent] = {
    labels.value map { l => MetadataEvent(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Labels}:${l.key}"), MetadataValue(l.value)) }
  }
}

class LabelManagerActor(serviceRegistryActor: ActorRef) extends Actor with ActorLogging with DefaultJsonProtocol {

  implicit val ec = context.dispatcher

  import WorkflowJsonSupport._

  def receive = {
    case LabelAddition(data) =>
      serviceRegistryActor ! PutMetadataActionAndRespond(labelsToMetadataEvents(data.labels, data.workflowId), self)
    case MetadataWriteSuccess(events) =>
      val wfId = events.head.key.workflowId
      val response = processLabelsResponse(wfId, metadataEventsToLabels(events))
      context.parent ! RequestComplete((StatusCodes.OK, response))
    case MetadataWriteFailure(failure, events) =>
      val wfId = events.head.key.workflowId
      val response = APIResponse.fail(new RuntimeException(s"Can't find metadata service to update labels for $wfId due to ${failure.getMessage}"))
      context.parent ! RequestComplete((StatusCodes.InternalServerError, response))
  }
}
