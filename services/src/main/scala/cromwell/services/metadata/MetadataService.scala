package cromwell.services.metadata

import java.time.OffsetDateTime

import io.circe.Json
import akka.actor.ActorRef
import cats.data.NonEmptyList
import cromwell.core._
import cromwell.services.ServiceRegistryActor.{ListenToMessage, ServiceRegistryMessage}
import common.exception.{MessageAggregation, ThrowableAggregation}
import wom.core._
import wom.values._

import scala.util.Random

object MetadataService {

  final val MetadataServiceName = "MetadataService"

  final case class WorkflowQueryResult(id: String,
                                       name: Option[String],
                                       status: Option[String],
                                       submission: Option[OffsetDateTime],
                                       start: Option[OffsetDateTime],
                                       end: Option[OffsetDateTime],
                                       labels: Option[Map[String, String]],
                                       parentWorkflowId: Option[String],
                                       rootWorkflowId: Option[String])

  final case class WorkflowQueryResponse(results: Seq[WorkflowQueryResult], totalResultsCount: Int)

  final case class QueryMetadata(page: Option[Int], pageSize: Option[Int], totalRecords: Option[Int])

  trait MetadataServiceMessage
  /**
    * Command Actions
    */
  trait MetadataServiceAction extends MetadataServiceMessage with ServiceRegistryMessage {
    def serviceName = MetadataServiceName
  }
  trait MetadataReadAction extends MetadataServiceAction

  trait WorkflowMetadataReadAction extends MetadataReadAction {
    def workflowId: WorkflowId
  }

  object PutMetadataAction {
    def apply(event: MetadataEvent, others: MetadataEvent*) = new PutMetadataAction(List(event) ++ others)
  }

  /**
    * Import from here with care! We extend every ActorRef, so import as locally as possible!
    */
  object implicits {
    implicit class MetadataAutoPutter(serviceRegistryActor: ActorRef) {
      def putMetadata(workflowId: WorkflowId, jobKey: Option[JobKey], keyValue: Map[String, Any]) = {
        val metadataJobKey = jobKey map { jk => MetadataJobKey(jk.node.fullyQualifiedName, jk.index, jk.attempt) }

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

  sealed trait MetadataWriteAction extends MetadataServiceAction {
    def maxAttempts: Int
    def events: Iterable[MetadataEvent]
    def size: Int = events.size
  }

  val MaximumMetadataWriteAttempts = 10
  final case class PutMetadataAction(events: Iterable[MetadataEvent], maxAttempts: Int = MaximumMetadataWriteAttempts) extends MetadataWriteAction
  final case class PutMetadataActionAndRespond(events: Iterable[MetadataEvent], replyTo: ActorRef, maxAttempts: Int = MaximumMetadataWriteAttempts) extends MetadataWriteAction

  final case object ListenToMetadataWriteActor extends MetadataServiceAction with ListenToMessage

  // Utility object to get GetMetadataAction's for a workflow-only query:
  object GetSingleWorkflowMetadataAction {
    def apply(workflowId: WorkflowId,
              includeKeysOption: Option[NonEmptyList[String]],
              excludeKeysOption: Option[NonEmptyList[String]],
              expandSubWorkflows: Boolean): WorkflowMetadataReadAction = {
      GetMetadataAction(MetadataQuery(workflowId, None, None, includeKeysOption, excludeKeysOption, expandSubWorkflows))
    }
  }


  final case class GetMetadataAction(key: MetadataQuery) extends WorkflowMetadataReadAction {
    override def workflowId: WorkflowId = key.workflowId
  }
  final case class GetStatus(workflowId: WorkflowId) extends WorkflowMetadataReadAction
  final case class GetLabels(workflowId: WorkflowId) extends WorkflowMetadataReadAction
  final case class QueryForWorkflowsMatchingParameters(parameters: Seq[(String, String)]) extends MetadataReadAction
  final case class WorkflowOutputs(workflowId: WorkflowId) extends WorkflowMetadataReadAction
  final case class GetLogs(workflowId: WorkflowId) extends WorkflowMetadataReadAction
  case object RefreshSummary extends MetadataServiceAction
  trait ValidationCallback {
    def onMalformed(possibleWorkflowId: String): Unit
    def onRecognized(workflowId: WorkflowId): Unit
    def onUnrecognized(possibleWorkflowId: String): Unit
    def onFailure(possibleWorkflowId: String, throwable: Throwable): Unit
  }

  final case class ValidateWorkflowIdInMetadata(possibleWorkflowId: WorkflowId) extends MetadataServiceAction
  final case class ValidateWorkflowIdInMetadataSummaries(possibleWorkflowId: WorkflowId) extends MetadataServiceAction

  /**
    * Responses
    */
  trait MetadataServiceResponse extends MetadataServiceMessage
  trait MetadataServiceFailure extends MetadataServiceResponse {
    def reason: Throwable
  }

  final case class MetadataLookupJsonResponse(query: MetadataQuery, result: Json) extends MetadataServiceResponse
  final case class MetadataLookupFailed(query: MetadataQuery, reason: Throwable)

  final case class MetadataLookupResponse(query: MetadataQuery, eventList: Seq[MetadataEvent]) extends MetadataServiceResponse
  final case class MetadataServiceKeyLookupFailed(query: MetadataQuery, reason: Throwable) extends MetadataServiceFailure

  final case class StatusLookupResponse(workflowId: WorkflowId, status: WorkflowState) extends MetadataServiceResponse
  final case class StatusLookupFailed(workflowId: WorkflowId, reason: Throwable) extends MetadataServiceFailure

  final case class LabelLookupResponse(workflowId: WorkflowId, labels: Map[String, String]) extends MetadataServiceResponse
  final case class LabelLookupFailed(workflowId: WorkflowId, reason: Throwable) extends MetadataServiceFailure

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

  private implicit class EnhancedWomTraversable(val womValues: Traversable[WomValue]) extends AnyVal {
    def toEvents(metadataKey: MetadataKey): List[MetadataEvent] = if (womValues.isEmpty) {
      List(MetadataEvent.empty(metadataKey.copy(key = s"${metadataKey.key}[]")))
    } else {
      womValues.toList
        .zipWithIndex
        .flatMap { case (value, index) => womValueToMetadataEvents(metadataKey.copy(key = s"${metadataKey.key}[$index]"), value) }
    }
  }
  
  private def toPrimitiveEvent(metadataKey: MetadataKey, valueName: String)(value: Option[Any]) = value match {
    case Some(v) => MetadataEvent(metadataKey.copy(key = s"${metadataKey.key}:$valueName"), MetadataValue(v))
    case None => MetadataEvent(metadataKey.copy(key = s"${metadataKey.key}:$valueName"), MetadataValue("", MetadataNull))
  }
  
  def womValueToMetadataEvents(metadataKey: MetadataKey, womValue: WomValue): Iterable[MetadataEvent] = womValue match {
    case WomArray(_, valueSeq) => valueSeq.toEvents(metadataKey)
    case WomMap(_, valueMap) =>
      if (valueMap.isEmpty) {
        List(MetadataEvent.empty(metadataKey))
      } else {
        valueMap.toList flatMap { case (key, value) => womValueToMetadataEvents(metadataKey.copy(key = metadataKey.key + s":${key.valueString}"), value) }
      }
    case objectLike: WomObjectLike =>
      if (objectLike.values.isEmpty) {
        List(MetadataEvent.empty(metadataKey))
      } else {
        objectLike.values.toList flatMap { case (key, value) => womValueToMetadataEvents(metadataKey.copy(key = metadataKey.key + s":$key"), value) }
      }
    case WomOptionalValue(_, Some(value)) =>
      womValueToMetadataEvents(metadataKey, value)
    case WomPair(left, right) =>
      womValueToMetadataEvents(metadataKey.copy(key = metadataKey.key + ":left"), left) ++
        womValueToMetadataEvents(metadataKey.copy(key = metadataKey.key + ":right"), right)
    case populated: WomMaybePopulatedFile =>
      import mouse.all._
      val secondaryFiles = populated.secondaryFiles.toEvents(metadataKey.copy(key = s"${metadataKey.key}:secondaryFiles"))

      List(
        MetadataEvent(metadataKey.copy(key = s"${metadataKey.key}:class"), MetadataValue("File")),
        populated.valueOption |> toPrimitiveEvent(metadataKey, "location"),
        populated.checksumOption |> toPrimitiveEvent(metadataKey, "checksum"),
        populated.sizeOption |> toPrimitiveEvent(metadataKey, "size"),
        populated.formatOption |> toPrimitiveEvent(metadataKey, "format"),
        populated.contentsOption |> toPrimitiveEvent(metadataKey, "contents")
      ) ++ secondaryFiles
    case listedDirectory: WomMaybeListedDirectory =>
      import mouse.all._
      val listing = listedDirectory.listingOption.toList.flatten.toEvents(metadataKey.copy(key = s"${metadataKey.key}:listing"))
      List(
        MetadataEvent(metadataKey.copy(key = s"${metadataKey.key}:class"), MetadataValue("Directory")),
        listedDirectory.valueOption |> toPrimitiveEvent(metadataKey, "location")
      ) ++ listing
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

      case _ =>
        val message = List(MetadataEvent(metadataKey.copy(key = s"$metadataKeyAndFailureIndex:message"), MetadataValue(t.getMessage)))
        val causeKey = metadataKey.copy(key = s"$metadataKeyAndFailureIndex:causedBy")
        val cause = Option(t.getCause) map { cause => throwableToMetadataEvents(causeKey, cause, 0) } getOrElse emptyCauseList
        message ++ cause
    }
  }
}

