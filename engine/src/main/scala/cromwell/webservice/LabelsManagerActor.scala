package cromwell.webservice

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core.labels.Labels
import cromwell.core._
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.services.metadata.MetadataService._
import cromwell.webservice.LabelsManagerActor._
import spray.json.{DefaultJsonProtocol, JsObject, JsString}

import scala.language.postfixOps


object LabelsManagerActor {

  def props(serviceRegistryActor: ActorRef) = Props(new LabelsManagerActor(serviceRegistryActor)).withDispatcher(Dispatcher.ApiDispatcher)

  final case class LabelsData(workflowId: WorkflowId, labels: Labels)

  sealed trait LabelsMessage {
    def data: LabelsData
  }

  sealed trait LabelsAction extends LabelsMessage
  final case class LabelsAddition(data: LabelsData) extends LabelsAction

  sealed trait LabelsResponse extends LabelsMessage

  def processLabelsResponse(workflowId: WorkflowId, labels: Map[String, String]): JsObject = {
    JsObject(Map(
      WorkflowMetadataKeys.Id -> JsString(workflowId.toString),
      WorkflowMetadataKeys.Labels -> JsObject(labels mapValues JsString.apply)
    ))
  }

  def metadataEventsToLabels(events: Iterable[MetadataEvent]): Map[String, String] = {
    events collect { case MetadataEvent(MetadataKey(_, _, key), Some(MetadataValue(value, _)), _) => key.split("\\:").last -> value } toMap
  }

  def labelsToMetadataEvents(labels: Labels, workflowId: WorkflowId): Iterable[MetadataEvent] = {
    labels.value map { l => MetadataEvent(MetadataKey(workflowId, None, s"${WorkflowMetadataKeys.Labels}:${l.key}"), MetadataValue(l.value)) }
  }

  sealed abstract class LabelsManagerActorResponse
  final case class BuiltLabelsManagerResponse(response: JsObject) extends LabelsManagerActorResponse
  final case class FailedLabelsManagerResponse(reason: Throwable) extends LabelsManagerActorResponse
}

class LabelsManagerActor(serviceRegistryActor: ActorRef) extends Actor with ActorLogging with DefaultJsonProtocol {

  implicit val ec = context.dispatcher

  private var wfId: Option[WorkflowId] = None
  private var target: ActorRef = ActorRef.noSender

  def receive = {
    case LabelsAddition(data) =>
      wfId = Option(data.workflowId)
      target = sender()
      serviceRegistryActor ! PutMetadataActionAndRespond(labelsToMetadataEvents(data.labels, data.workflowId), self)
    case MetadataWriteSuccess(events) =>
      target ! BuiltLabelsManagerResponse(processLabelsResponse(wfId.get, metadataEventsToLabels(events)))
      context stop self
    case MetadataWriteFailure(failure, events @ _) =>
      target ! FailedLabelsManagerResponse(new RuntimeException(s"Unable to update labels for ${wfId.get} due to ${failure.getMessage}"))
      context stop self
  }
}
