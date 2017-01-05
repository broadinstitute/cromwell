package cromwell.services.metadata

import java.time.OffsetDateTime

import akka.actor.{ActorRef, DeadLetterSuppression}
import cats.data.NonEmptyList
import cromwell.core.{JobKey, WorkflowId, WorkflowState}
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import lenthall.exception.ThrowableAggregation
import wdl4s.values._


object MetadataService {

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
  trait ReadAction extends MetadataServiceAction
  object PutMetadataAction {
    def apply(event: MetadataEvent, others: MetadataEvent*) = new PutMetadataAction(List(event) ++ others)
  }

  /**
    * Import from here with care! We extend every ActorRef, so import as locally as possible!
    */
  object implicits {
    implicit class MetadataAutoPutter(serviceRegistryActor: ActorRef) {
      def putMetadata(workflowId: WorkflowId, jobKey: Option[JobKey], keyValue: Map[String, Any]) = {
        val metadataJobKey = jobKey map { jk => MetadataJobKey(jk.scope.fullyQualifiedName, jk.index, jk.attempt) }

        val events = keyValue map { case (key, value) =>
          val metadataKey = MetadataKey(workflowId, metadataJobKey, key)
          MetadataEvent(metadataKey, MetadataValue(value))
        }
        serviceRegistryActor ! PutMetadataAction(events)
      }
    }
  }

  case class PutMetadataAction(events: Iterable[MetadataEvent]) extends MetadataServiceAction
  case class GetSingleWorkflowMetadataAction(workflowId: WorkflowId, includeKeysOption: Option[NonEmptyList[String]],
                                             excludeKeysOption: Option[NonEmptyList[String]],
                                             expandSubWorkflows: Boolean)
    extends ReadAction
  case class GetMetadataQueryAction(key: MetadataQuery) extends ReadAction
  case class GetStatus(workflowId: WorkflowId) extends ReadAction
  case class WorkflowQuery[A](uri: A, parameters: Seq[(String, String)]) extends ReadAction
  case class WorkflowOutputs(workflowId: WorkflowId) extends ReadAction
  case class GetLogs(workflowId: WorkflowId) extends ReadAction
  case object RefreshSummary extends MetadataServiceAction
  trait ValidationCallback {
    def onMalformed(possibleWorkflowId: String): Unit
    def onRecognized(workflowId: WorkflowId): Unit
    def onUnrecognized(possibleWorkflowId: String): Unit
    def onFailure(possibleWorkflowId: String, throwable: Throwable): Unit
  }
  final case class ValidateWorkflowIdAndExecute(possibleWorkflowId: String,
                                                validationCallback: ValidationCallback) extends MetadataServiceAction

  /**
    * Responses
    */
  trait MetadataServiceResponse extends MetadataServiceMessage
  trait MetadataServiceFailure extends MetadataServiceResponse {
    def reason: Throwable
  }

  case class MetadataPutAcknowledgement(putRequest: PutMetadataAction) extends MetadataServiceResponse with DeadLetterSuppression
  case class MetadataPutFailed(putRequest: PutMetadataAction, reason: Throwable) extends MetadataServiceFailure

  case class MetadataLookupResponse(query: MetadataQuery, eventList: Seq[MetadataEvent]) extends MetadataServiceResponse
  case class MetadataServiceKeyLookupFailed(query: MetadataQuery, reason: Throwable) extends MetadataServiceFailure

  case class StatusLookupResponse(workflowId: WorkflowId, status: WorkflowState) extends MetadataServiceResponse
  case class StatusLookupFailed(workflowId: WorkflowId, reason: Throwable) extends MetadataServiceFailure

  final case class WorkflowQuerySuccess[A](uri: A, response: WorkflowQueryResponse, meta: Option[QueryMetadata])
    extends MetadataServiceResponse
  final case class WorkflowQueryFailure(reason: Throwable) extends MetadataServiceFailure

  final case class WorkflowOutputsResponse(id: WorkflowId, outputs: Seq[MetadataEvent]) extends MetadataServiceResponse
  final case class WorkflowOutputsFailure(id: WorkflowId, reason: Throwable) extends MetadataServiceFailure

  final case class LogsResponse(id: WorkflowId, logs: Seq[MetadataEvent]) extends MetadataServiceResponse
  final case class LogsFailure(id: WorkflowId, reason: Throwable) extends MetadataServiceFailure

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
    case WdlOptionalValue(_, Some(value)) =>
      wdlValueToMetadataEvents(metadataKey, value)
    case WdlPair(left, right) =>
      wdlValueToMetadataEvents(metadataKey.copy(key = metadataKey.key + ":left"), left) ++
        wdlValueToMetadataEvents(metadataKey.copy(key = metadataKey.key + ":right"), right)
    case value =>
      List(MetadataEvent(metadataKey, MetadataValue(value)))
  }

  def throwableToMetadataEvents(metadataKey: MetadataKey, t: Throwable): List[MetadataEvent] = {
    t match {
      case aggregation: ThrowableAggregation =>
        aggregation.errorMessages.toList map { message => MetadataEvent(metadataKey.copy(key = s"${metadataKey.key}:message"), MetadataValue(s"${aggregation.exceptionContext}: $message")) }
      case other =>
        val message = List(MetadataEvent(metadataKey.copy(key = s"${metadataKey.key}:message"), MetadataValue(t.getMessage)))
        val cause = Option(t.getCause) map { cause => throwableToMetadataEvents(metadataKey.copy(key = s"${metadataKey.key}:causedBy"), cause) } getOrElse List.empty
        message ++ cause
    }
  }
}

