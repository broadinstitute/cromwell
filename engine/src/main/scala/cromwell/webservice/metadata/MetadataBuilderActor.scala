package cromwell.webservice.metadata

import java.util.UUID

import akka.actor.{ActorRef, LoggingFSM, Props}
import cats.effect.IO
import cromwell.core.Dispatcher.ApiDispatcher
import cromwell.core._
import cromwell.services.ServiceRegistryActor.ServiceRegistryFailure
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.webservice.metadata.MetadataBuilderActor._
import org.slf4j.LoggerFactory
import spray.json._


object MetadataBuilderActor {
  sealed abstract class MetadataBuilderActorResponse
  case class BuiltMetadataResponse(response: JsObject) extends MetadataBuilderActorResponse
  case class FailedMetadataResponse(reason: Throwable) extends MetadataBuilderActorResponse

  sealed trait MetadataBuilderActorState
  case object Idle extends MetadataBuilderActorState
  case object WaitingForMetadataService extends MetadataBuilderActorState
  case object WaitingForSubWorkflows extends MetadataBuilderActorState

  case class MetadataBuilderActorData(
                                       originalQuery: MetadataQuery,
                                       originalEvents: Seq[MetadataEvent],
                                       subWorkflowsMetadata: Map[String, JsValue],
                                       waitFor: Int
                                     ) {
    def withSubWorkflow(id: String, metadata: JsValue) = {
      this.copy(subWorkflowsMetadata = subWorkflowsMetadata + ((id, metadata)))
    }

    def isComplete = subWorkflowsMetadata.size == waitFor
  }

  def props(serviceRegistryActor: ActorRef) = {
    Props(new MetadataBuilderActor(serviceRegistryActor)).withDispatcher(ApiDispatcher)
  }

  val log = LoggerFactory.getLogger("MetadataBuilder")

  def uniqueActorName: String = List("MetadataBuilderActor", UUID.randomUUID()).mkString("-")
}

class MetadataBuilderActor(serviceRegistryActor: ActorRef) extends LoggingFSM[MetadataBuilderActorState, Unit]
  with DefaultJsonProtocol {
  import MetadataBuilderActor._

  private var target: ActorRef = ActorRef.noSender

  implicit val ec = context.dispatcher
  
  startWith(Idle, ())
  val tag = self.path.name

  when(Idle) {
    case Event(GetSingleWorkflowMetadataAction(workflowId, includeKeysOption, excludeKeysOption, expandSubWorkflows), _) =>
      target = sender()
      val includeKeys = if (expandSubWorkflows)
        includeKeysOption map { _.::(CallMetadataKeys.SubWorkflowId) }
      else 
        includeKeysOption
      buildAndStop(MetadataQuery(workflowId, None, None, includeKeys, excludeKeysOption, expandSubWorkflows))
    case Event(GetMetadataQueryAction(query@MetadataQuery(_, _, _, _, _, _)), _) =>
      target = sender()
      buildAndStop(query)
    case Event(GetLogs(workflowId), _) =>
      target = sender()
      buildAndStop(StreamMetadataBuilder.workflowLogs(workflowId))
    case Event(WorkflowOutputs(workflowId), _) =>
      target = sender()
      buildAndStop(StreamMetadataBuilder.workflowOutputs(workflowId))
    case Event(action: MetadataServiceAction, _) =>
      target = sender()
      serviceRegistryActor ! action
      goto(WaitingForMetadataService)
  }

  private def allDone = {
    context stop self
    stay()
  }

  when(WaitingForMetadataService) {
    case Event(StatusLookupResponse(w, status), _) =>
      target ! BuiltMetadataResponse(processStatusResponse(w, status))
      allDone
    case Event(LabelLookupResponse(w, labels), _) =>
      target ! BuiltMetadataResponse(processLabelsResponse(w, labels))
      allDone
    case Event(_: ServiceRegistryFailure, _) =>
      target ! FailedMetadataResponse(new RuntimeException("Can't find metadata service"))
      allDone
    case Event(failure: MetadataServiceFailure, _) =>
      target ! FailedMetadataResponse(failure.reason)
      allDone
    case Event(unexpectedMessage, stateData) =>
      target ! FailedMetadataResponse(new RuntimeException(s"MetadataBuilderActor $tag(WaitingForMetadataService, $stateData) got an unexpected message: $unexpectedMessage"))
      context stop self
      stay()
  }

  whenUnhandled {
    case Event(message, data) =>
      log.error(s"Received unexpected message $message in state $stateName with data $data")
      stay()
  }

  def failAndDie(reason: Throwable) = {
    target ! FailedMetadataResponse(reason)
    context stop self
    stay()
  }

  def buildAndStop(ioComponent: IO[MetadataComponent]): State = {
    target ! BuiltMetadataResponse(ioComponent.unsafeRunSync().toJson.asJsObject)
    log.info("Done")
    allDone
  }

  def buildAndStop(query: MetadataQuery): State = {
    target ! BuiltMetadataResponse(StreamMetadataBuilder.workflowMetadataQuery(query).unsafeRunSync())
    log.info("Done")
    allDone
  }

  def processStatusResponse(workflowId: WorkflowId, status: WorkflowState): JsObject = {
    JsObject(Map(
      WorkflowMetadataKeys.Status -> JsString(status.toString),
      WorkflowMetadataKeys.Id -> JsString(workflowId.toString)
    ))
  }

  def processLabelsResponse(workflowId: WorkflowId, labels: Map[String, String]): JsObject = {
    val jsLabels = labels map { case (k, v) => k -> JsString(v) }
    JsObject(Map(
      WorkflowMetadataKeys.Id -> JsString(workflowId.toString),
      WorkflowMetadataKeys.Labels -> JsObject(jsLabels)
    ))
  }
}
