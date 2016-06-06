package cromwell.services

import java.time.OffsetDateTime

import akka.actor.ActorRef
import cromwell.core.{WorkflowId, WorkflowState}
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import spray.http.Uri
import wdl4s.values._

import scala.language.postfixOps

object MetadataServiceActor {

  final val MetadataServiceName = "MetadataService"

  case class WorkflowQueryResult(id: String, name: Option[String], status: Option[String], start: Option[OffsetDateTime], end: Option[OffsetDateTime])

  case class WorkflowQueryResponse(results: Seq[WorkflowQueryResult])

  case class QueryMetadata(page: Option[Int], pageSize: Option[Int], totalRecords: Option[Int])

  trait MetadataServiceMessage
  /**
    * Command Actions
    */
  trait MetadataServiceAction extends MetadataServiceMessage with ServiceRegistryMessage {
    def serviceName = MetadataServiceName
  }
  object PutMetadataAction {
    def apply(event: MetadataEvent, others: MetadataEvent*) = new PutMetadataAction(List(event) ++ others)
  }
  case class PutMetadataAction(events: Iterable[MetadataEvent]) extends MetadataServiceAction
  case class GetSingleWorkflowMetadataAction(workflowId: WorkflowId) extends MetadataServiceAction
  case class GetMetadataQueryAction(key: MetadataQuery) extends MetadataServiceAction
  case class GetStatus(workflowId: WorkflowId) extends MetadataServiceAction
  case class WorkflowQuery(uri: Uri, parameters: Seq[(String, String)]) extends MetadataServiceAction
  case object RefreshSummary extends MetadataServiceAction
  final case class HandleNotFound(workflowId: WorkflowId, sndr: ActorRef) extends MetadataServiceAction

  /**
    * Responses
    */
  trait MetadataServiceResponse extends MetadataServiceMessage
  case class MetadataPutAcknowledgement(putRequest: PutMetadataAction) extends MetadataServiceResponse
  case class MetadataPutFailed(putRequest: PutMetadataAction, reason: Throwable) extends MetadataServiceResponse

  case class MetadataLookupResponse(query: MetadataQuery, eventList: Seq[MetadataEvent]) extends MetadataServiceResponse
  case class MetadataServiceKeyLookupFailed(query: MetadataQuery, reason: Throwable) extends MetadataServiceResponse

  case class StatusLookupResponse(workflowId: WorkflowId, status: WorkflowState) extends MetadataServiceResponse
  case class StatusLookupNotFound(workflowId: WorkflowId) extends MetadataServiceResponse
  case class StatusLookupFailed(workflowId: WorkflowId, reason: Throwable) extends MetadataServiceResponse

  final case class WorkflowQuerySuccess(uri: Uri, response: WorkflowQueryResponse, meta: Option[QueryMetadata]) extends MetadataServiceResponse
  final case class WorkflowQueryFailure(failure: Throwable) extends MetadataServiceResponse

  /* TODO: PBE: No EngineMetadataServiceActor.props until circular dependencies fixed.
  def props(serviceConfig: Config, globalConfig: Config) = {
    Props(MetadataServiceActor(serviceConfig, globalConfig))
  }
  */

  def wdlValueToMetadataEvents(metadataKey: MetadataKey, wdlValue: WdlValue): Iterable[MetadataEvent] = wdlValue match {
    case WdlArray(_, valueSeq) =>
      if (valueSeq.isEmpty) {
        List(MetadataEvent.empty(metadataKey.copy(key = s"${metadataKey.key}[]")))
      } else {
        val zippedSeq = valueSeq.zipWithIndex
        zippedSeq.toList flatMap { case (value, index) => wdlValueToMetadataEvents(metadataKey.copy(key = s"${metadataKey.key}[$index]"), value) }
      }
    case WdlMap(_, valueMap) =>
      if (valueMap.isEmpty) {
        List(MetadataEvent.empty(metadataKey))
      } else {
        valueMap.toList flatMap { case (key, value) => wdlValueToMetadataEvents(metadataKey.copy(key = metadataKey.key + s":${key.valueString}"), value) }
      }
    case value =>
      List(MetadataEvent(metadataKey, MetadataValue(value)))
  }

  def throwableToMetadataEvents(metadataKey: MetadataKey, t: Throwable): List[MetadataEvent] = {
    val message = List(MetadataEvent(metadataKey.copy(key = s"${metadataKey.key}:message"), MetadataValue(t.getMessage)))
    val cause = Option(t.getCause) map { cause => throwableToMetadataEvents(metadataKey.copy(key = s"${metadataKey.key}:causedBy"), cause) } getOrElse List.empty
    message ++ cause
  }
}

object MetadataServiceActorImplicits {
  import MetadataServiceActor._

  implicit class EnhancedServiceRegistryActorForMetadata(val actor: ActorRef) extends AnyVal {
    def pushWdlValueMetadata(metadataKey: MetadataKey, output: WdlValue): Unit = actor ! PutMetadataAction(wdlValueToMetadataEvents(metadataKey, output))
  }
}
