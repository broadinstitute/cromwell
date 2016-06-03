package cromwell.services

import java.time.OffsetDateTime

import akka.actor.ActorRef
import cromwell.core.{WorkflowId, WorkflowState}
import cromwell.services.MetadataServiceActor.PutMetadataAction
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
  case class PutMetadataAction(event: MetadataEvent) extends MetadataServiceAction
  case class GetSingleWorkflowMetadataAction(workflowId: WorkflowId) extends MetadataServiceAction
  case class GetMetadataQueryAction(key: MetadataQuery) extends MetadataServiceAction
  case class GetStatus(workflowId: WorkflowId) extends MetadataServiceAction
  case class WorkflowQuery(uri: Uri, parameters: Seq[(String, String)]) extends MetadataServiceAction
  case object RefreshSummary extends MetadataServiceAction

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
}

object MetadataServiceActorImplicits {
  implicit class EnhancedServiceRegistryActorForMetadata(val actor: ActorRef) extends AnyVal {
    def pushWdlValueMetadata(metadataKey: MetadataKey, output: WdlValue): Unit = output match {
      case WdlArray(_, valueSeq) =>
        val zippedSeq = valueSeq.zipWithIndex
        zippedSeq foreach { case (value, index) => actor.pushWdlValueMetadata(metadataKey.copy(key = s"${metadataKey.key}[$index]"), value) }
      case WdlMap(_, valueMap) =>
        valueMap foreach { case (key, value) => actor.pushWdlValueMetadata(metadataKey.copy(key = metadataKey.key + s":${key.valueString}"), value) }
      case value =>
        actor ! PutMetadataAction(MetadataEvent(metadataKey, MetadataValue(value)))
    }

    def pushThrowableMetadata(metadataKey: MetadataKey, t: Throwable): Unit = {
      actor ! PutMetadataAction(MetadataEvent(metadataKey.copy(key = s"${metadataKey.key}:message"), MetadataValue(t.getMessage)))
      Option(t.getCause) foreach { cause => pushThrowableMetadata(metadataKey.copy(key = s"${metadataKey.key}:causedBy"), cause) }
    }
  }
}
