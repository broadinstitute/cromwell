package cromwell.webservice

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import cromwell.core._
import cromwell.core.labels.{Label, Labels}
import cromwell.services.metadata.MetadataEvent
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.BuiltMetadataResponse
import cromwell.webservice.LabelsManagerActor._
import spray.json.{DefaultJsonProtocol, JsObject, JsString}

object LabelsManagerActor {

  def props(serviceRegistryActor: ActorRef): Props = Props(new LabelsManagerActor(serviceRegistryActor)).withDispatcher(Dispatcher.ApiDispatcher)

  final case class LabelsData(workflowId: WorkflowId, labels: Labels)

  sealed trait LabelsMessage {
    def data: LabelsData
  }

  sealed trait LabelsAction extends LabelsMessage
  final case class LabelsAddition(data: LabelsData) extends LabelsAction

  sealed trait LabelsResponse extends LabelsMessage

  sealed abstract class LabelsManagerActorResponse
  final case class BuiltLabelsManagerResponse(response: JsObject) extends LabelsManagerActorResponse
  final case class FailedLabelsManagerResponse(reason: Throwable) extends LabelsManagerActorResponse
}

class LabelsManagerActor(serviceRegistryActor: ActorRef) extends Actor with ActorLogging with DefaultJsonProtocol {

  implicit val ec = context.dispatcher

  // These var Options get set by the effective entry point of the actor, the LabelsAddition case
  private var wfId: Option[WorkflowId] = None
  private var newLabels: Option[Labels] = None

  private var target: ActorRef = ActorRef.noSender

  override def receive = {
    case LabelsAddition(data) =>
      wfId = Option(data.workflowId)
      newLabels = Option(data.labels)
      target = sender()
      val metadataEvents = MetadataEvent.labelsToMetadataEvents(data.labels, data.workflowId)
      serviceRegistryActor ! PutMetadataActionAndRespond(metadataEvents, self)
    case MetadataWriteSuccess(_) =>
      /*
        Ask the metadata store for the current set of labels, so we can return the full label set to the user.
        At this point in the actor lifecycle, wfId has already been filled out so the .get is safe
      */
      serviceRegistryActor ! GetLabels(wfId.get)
    case BuiltMetadataResponse(_, jsObject) =>
      /*
        There's some trickery going on here. We've updated the labels in the metadata store but almost certainly when
        the store received the GetLabels request above the summarizer will not have been run so our new values are
        not present. Instead we'll fake it by manually merging what the metadata store has, as well as the new updated
        values.

        There's a potential flaw here in that a second PATCH request to the labels endpoint could come in prior to
        summarization and thus the response to the user would lose the previous updates, but we can cross that bridge
        if it becomes a problem. So if you find yourself debugging why previously updated labels aren't showing up in
        the return packet, this is a likely cause.

        At this point in the actor lifecycle, newLabels will have been filled in so the .get is safe
      */

      def replaceOrAddLabel(originalJson: JsObject, label: Label): JsObject = {
        val labels = originalJson.fields.get("labels").map(_.asJsObject.fields).getOrElse(Map.empty)
        val updatedLabels = labels + (label.key -> JsString(label.value))

        JsObject(originalJson.fields + ("labels" -> JsObject(updatedLabels)))
      }

      val updatedJson = newLabels.get.value.foldLeft(jsObject)(replaceOrAddLabel)
      target ! BuiltLabelsManagerResponse(updatedJson)
      context stop self
    case f: MetadataServiceFailure =>
      /*
        This case handles both a MetadataWriteFailure from the initial write as well as a LabelLookupFailure from
        the label lookup for the user response.

        At this point in the actor lifecycle, wfId has already been filled out so the .get is safe
       */
      target ! FailedLabelsManagerResponse(new RuntimeException(s"Unable to update labels for ${wfId.get} due to ${f.reason.getMessage}"))
      context stop self
  }
}
