package cromwell.services.metadata

import java.time.OffsetDateTime

import akka.actor.ActorRef
import cats.data.NonEmptyList
import cromwell.core.{FullyQualifiedName, JobKey, WorkflowId, WorkflowState}
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import lenthall.exception.{MessageAggregation, ThrowableAggregation}
import wdl4s.values._

import scala.util.Random


object MetadataService {

  final val MetadataServiceName = "MetadataService"

  final case class WorkflowQueryResult(id: String, name: Option[String], status: Option[String], start: Option[OffsetDateTime], end: Option[OffsetDateTime])

  final case class WorkflowQueryResponse(results: Seq[WorkflowQueryResult])

  final case class QueryMetadata(page: Option[Int], pageSize: Option[Int], totalRecords: Option[Int])

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

      def putMetadataWithRawKey(workflowId: WorkflowId, jobKey: Option[(FullyQualifiedName, Option[Int], Int)], keyValue: Map[String, Any]) = {
        val metadataJobKey = jobKey map { case (fullyQualifiedName, index, attempt) => MetadataJobKey(fullyQualifiedName, index, attempt) }

        val events = keyValue map { case (key, value) =>
          val metadataKey = MetadataKey(workflowId, metadataJobKey, key)
          MetadataEvent(metadataKey, MetadataValue(value))
        }
        serviceRegistryActor ! PutMetadataAction(events)
      }
    }
  }

  final case class PutMetadataAction(events: Iterable[MetadataEvent]) extends MetadataServiceAction
  final case class PutMetadataActionAndRespond(events: Iterable[MetadataEvent], replyTo: ActorRef) extends MetadataServiceAction

  case class GetSingleWorkflowMetadataAction(workflowId: WorkflowId, includeKeysOption: Option[NonEmptyList[String]],
                                             excludeKeysOption: Option[NonEmptyList[String]],
                                             expandSubWorkflows: Boolean)
    extends ReadAction
  case class GetMetadataQueryAction(key: MetadataQuery) extends ReadAction
  case class GetStatus(workflowId: WorkflowId) extends ReadAction
  case class WorkflowQuery(parameters: Seq[(String, String)]) extends ReadAction
  case class WorkflowOutputs(workflowId: WorkflowId) extends ReadAction
  case class GetLogs(workflowId: WorkflowId) extends ReadAction
  case object RefreshSummary extends MetadataServiceAction
  trait ValidationCallback {
    def onMalformed(possibleWorkflowId: String): Unit
    def onRecognized(workflowId: WorkflowId): Unit
    def onUnrecognized(possibleWorkflowId: String): Unit
    def onFailure(possibleWorkflowId: String, throwable: Throwable): Unit
  }

  final case class ValidateWorkflowId(possibleWorkflowId: WorkflowId) extends MetadataServiceAction

  /**
    * Responses
    */
  trait MetadataServiceResponse extends MetadataServiceMessage
  trait MetadataServiceFailure extends MetadataServiceResponse {
    def reason: Throwable
  }

  final case class MetadataLookupResponse(query: MetadataQuery, eventList: Seq[MetadataEvent]) extends MetadataServiceResponse
  final case class MetadataServiceKeyLookupFailed(query: MetadataQuery, reason: Throwable) extends MetadataServiceFailure

  final case class StatusLookupResponse(workflowId: WorkflowId, status: WorkflowState) extends MetadataServiceResponse
  final case class StatusLookupFailed(workflowId: WorkflowId, reason: Throwable) extends MetadataServiceFailure

  final case class WorkflowOutputsResponse(id: WorkflowId, outputs: Seq[MetadataEvent]) extends MetadataServiceResponse
  final case class WorkflowOutputsFailure(id: WorkflowId, reason: Throwable) extends MetadataServiceFailure

  final case class LogsResponse(id: WorkflowId, logs: Seq[MetadataEvent]) extends MetadataServiceResponse
  final case class LogsFailure(id: WorkflowId, reason: Throwable) extends MetadataServiceFailure

  final case class MetadataWriteSuccess(events: Iterable[MetadataEvent]) extends MetadataServiceResponse
  final case class MetadataWriteFailure(reason: Throwable, events: Iterable[MetadataEvent]) extends MetadataServiceFailure

  sealed abstract class WorkflowValidationResponse extends MetadataServiceResponse
  case object RecognizedWorkflowId extends WorkflowValidationResponse
  case object UnrecognizedWorkflowId extends WorkflowValidationResponse
  final case class FailedToCheckWorkflowId(cause: Throwable) extends WorkflowValidationResponse

  sealed abstract class MetadataQueryResponse extends MetadataServiceResponse
  final case class WorkflowQuerySuccess(response: WorkflowQueryResponse, meta: Option[QueryMetadata]) extends MetadataQueryResponse
  final case class WorkflowQueryFailure(reason: Throwable) extends MetadataQueryResponse

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

  def throwableToMetadataEvents(metadataKey: MetadataKey, t: Throwable, failureIndex: Int = Random.nextInt(Int.MaxValue)): List[MetadataEvent] = {
    val emptyCauseList = List(MetadataEvent.empty(metadataKey.copy(key = metadataKey.key + s"[$failureIndex]:causedBy[]")))
    val metadataKeyAndFailureIndex = s"${metadataKey.key}[$failureIndex]"

    t match {
      case aggregation: ThrowableAggregation =>
        val message = List(MetadataEvent(metadataKey.copy(key = s"$metadataKeyAndFailureIndex:message"), MetadataValue(aggregation.exceptionContext)))
        val indexedCauses = aggregation.throwables.toList.zipWithIndex
        val indexedCauseEvents = if (indexedCauses.nonEmpty) {
          indexedCauses flatMap { case (cause, index) =>
            val causeKey = metadataKey.copy(key = s"$metadataKeyAndFailureIndex:causedBy")
            throwableToMetadataEvents(causeKey, cause, index)
          }
        } else {
          emptyCauseList
        }
        message ++ indexedCauseEvents
      case aggregation: MessageAggregation =>
        val message = List(MetadataEvent(metadataKey.copy(key = s"$metadataKeyAndFailureIndex:message"), MetadataValue(aggregation.exceptionContext)))
        val indexedCauses = aggregation.errorMessages.toList.zipWithIndex
        val indexedCauseEvents = if (indexedCauses.nonEmpty) {
          indexedCauses flatMap { case (cause, index) =>
            val causeMessageKey = metadataKey.copy(key = s"$metadataKeyAndFailureIndex:causedBy[$index]:message")
            val causeCausedByKey = metadataKey.copy(key = s"$metadataKeyAndFailureIndex:causedBy[$index]:causedBy[]")
            List(MetadataEvent(causeMessageKey, MetadataValue(cause)), MetadataEvent.empty(causeCausedByKey))
          }
        } else {
          emptyCauseList
        }
        message ++ indexedCauseEvents

      case other =>
        val message = List(MetadataEvent(metadataKey.copy(key = s"$metadataKeyAndFailureIndex:message"), MetadataValue(t.getMessage)))
        val causeKey = metadataKey.copy(key = s"$metadataKeyAndFailureIndex:causedBy")
        val cause = Option(t.getCause) map { cause => throwableToMetadataEvents(causeKey, cause, 0) } getOrElse emptyCauseList
        message ++ cause
    }
  }
}

