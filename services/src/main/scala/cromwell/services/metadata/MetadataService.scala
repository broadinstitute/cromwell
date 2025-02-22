package cromwell.services.metadata

import java.time.OffsetDateTime
import akka.actor.ActorRef
import cats.data.NonEmptyList
import cromwell.core._
import cromwell.services.ServiceRegistryActor.{ListenToMessage, ServiceRegistryMessage}
import common.exception.{MessageAggregation, ThrowableAggregation}
import cromwell.core.path.Path
import cromwell.database.sql.tables.MetadataEntry
import slick.basic.DatabasePublisher
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
                                       rootWorkflowId: Option[String],
                                       metadataArchiveStatus: MetadataArchiveStatus
  )

  final case class WorkflowQueryResponse(results: Seq[WorkflowQueryResult], totalResultsCount: Int)

  final case class QueryMetadata(page: Option[Int], pageSize: Option[Int], totalRecords: Option[Int])

  trait MetadataServiceMessage

  /**
    * Command Actions
    */
  trait MetadataServiceAction extends MetadataServiceMessage with ServiceRegistryMessage {
    def serviceName = MetadataServiceName
  }
  sealed trait BuildMetadataJsonAction extends MetadataServiceAction

  sealed trait BuildWorkflowMetadataJsonAction extends BuildMetadataJsonAction {
    def workflowId: WorkflowId
  }
  sealed trait BuildWorkflowMetadataJsonWithOverridableSourceAction extends BuildWorkflowMetadataJsonAction

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

      def putMetadataWithRawKey(workflowId: WorkflowId,
                                jobKey: Option[(FullyQualifiedName, Option[Int], Int)],
                                keyValue: Map[String, Any]
      ) = {
        val metadataJobKey = jobKey map { case (fullyQualifiedName, index, attempt) =>
          MetadataJobKey(fullyQualifiedName, index, attempt)
        }

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

  val MaximumMetadataActionAttempts = 10
  final case class PutMetadataAction(events: Iterable[MetadataEvent], maxAttempts: Int = MaximumMetadataActionAttempts)
      extends MetadataWriteAction
  final case class PutMetadataActionAndRespond(events: Iterable[MetadataEvent],
                                               replyTo: ActorRef,
                                               maxAttempts: Int = MaximumMetadataActionAttempts
  ) extends MetadataWriteAction

  final case object ListenToMetadataWriteActor extends MetadataServiceAction with ListenToMessage

  // Utility object to get GetMetadataAction's for a workflow-only query:
  object GetSingleWorkflowMetadataAction {
    def apply(workflowId: WorkflowId,
              includeKeysOption: Option[NonEmptyList[String]],
              excludeKeysOption: Option[NonEmptyList[String]],
              expandSubWorkflows: Boolean
    ): BuildWorkflowMetadataJsonAction =
      GetMetadataAction(MetadataQuery(workflowId, None, None, includeKeysOption, excludeKeysOption, expandSubWorkflows))
  }

  final case class GetMetadataAction(key: MetadataQuery, checkTotalMetadataRowNumberBeforeQuerying: Boolean = true)
      extends BuildWorkflowMetadataJsonWithOverridableSourceAction {

    override def workflowId: WorkflowId = key.workflowId
  }

  final case class GetMetadataStreamAction(workflowId: WorkflowId) extends MetadataServiceAction

  final case class GetStatus(workflowId: WorkflowId) extends BuildWorkflowMetadataJsonAction
  final case class GetLabels(workflowId: WorkflowId) extends BuildWorkflowMetadataJsonAction
  final case class GetRootAndSubworkflowLabels(workflowId: WorkflowId) extends BuildWorkflowMetadataJsonAction
  final case class QueryForWorkflowsMatchingParameters(parameters: Seq[(String, String)])
      extends BuildMetadataJsonAction
  final case class WorkflowOutputs(workflowId: WorkflowId) extends BuildWorkflowMetadataJsonWithOverridableSourceAction
  final case class GetLogs(workflowId: WorkflowId) extends BuildWorkflowMetadataJsonWithOverridableSourceAction
  final case class GetCost(workflowId: WorkflowId) extends BuildWorkflowMetadataJsonWithOverridableSourceAction
  case object RefreshSummary extends MetadataServiceAction
  case object SendMetadataTableSizeMetrics extends MetadataServiceAction

  final case class ValidateWorkflowIdInMetadataSummaries(possibleWorkflowId: WorkflowId) extends MetadataServiceAction
  final case class FetchWorkflowMetadataArchiveStatusAndEndTime(workflowId: WorkflowId) extends MetadataServiceAction
  final case class FetchFailedJobsMetadataWithWorkflowId(workflowId: WorkflowId) extends BuildWorkflowMetadataJsonAction

  /**
    * Responses
    */
  trait MetadataServiceResponse extends MetadataServiceMessage
  trait MetadataServiceFailure extends MetadataServiceResponse {
    def reason: Throwable
  }

  final case class MetadataLookupStreamSuccess(id: WorkflowId, result: DatabasePublisher[MetadataEntry])
      extends MetadataServiceResponse
  final case class MetadataLookupStreamFailed(id: WorkflowId, reason: Throwable) extends MetadataServiceResponse
  final case class MetadataLookupFailedTooLargeResponse(query: MetadataQuery, metadataSizeRows: Int)
      extends MetadataServiceResponse
  final case class MetadataLookupFailedTimeoutResponse(query: MetadataQuery) extends MetadataServiceResponse
  final case class FetchFailedTasksTimeoutResponse(workflowId: WorkflowId) extends MetadataServiceResponse
  final case class MetadataLookupResponse(query: MetadataQuery, eventList: Seq[MetadataEvent])
      extends MetadataServiceResponse
  final case class FetchFailedJobsMetadataLookupResponse(events: Seq[MetadataEvent]) extends MetadataServiceResponse
  final case class FetchFailedJobsMetadataLookupFailed(workflowId: WorkflowId, reason: Throwable)
      extends MetadataServiceFailure
  final case class MetadataServiceKeyLookupFailed(query: MetadataQuery, reason: Throwable)
      extends MetadataServiceFailure

  final case class StatusLookupResponse(workflowId: WorkflowId, status: WorkflowState) extends MetadataServiceResponse
  final case class StatusLookupFailed(workflowId: WorkflowId, reason: Throwable) extends MetadataServiceFailure

  final case class LabelLookupResponse(workflowId: WorkflowId, labels: Map[String, String])
      extends MetadataServiceResponse
  final case class LabelLookupFailed(workflowId: WorkflowId, reason: Throwable) extends MetadataServiceFailure

  final case class RootAndSubworkflowLabelsLookupResponse(rootWorkflowId: WorkflowId,
                                                          labels: Map[WorkflowId, Map[String, String]]
  ) extends MetadataServiceResponse
  final case class RootAndSubworkflowLabelsLookupFailed(rootWorkflowId: WorkflowId, reason: Throwable)
      extends MetadataServiceFailure

  final case class WorkflowOutputsResponse(id: WorkflowId, outputs: Seq[MetadataEvent]) extends MetadataServiceResponse
  final case class WorkflowOutputsFailure(id: WorkflowId, reason: Throwable) extends MetadataServiceFailure

  final case class LogsResponse(id: WorkflowId, logs: Seq[MetadataEvent]) extends MetadataServiceResponse
  final case class LogsFailure(id: WorkflowId, reason: Throwable) extends MetadataServiceFailure

  final case class CostResponse(id: WorkflowId, status: WorkflowState, metadataResponse: MetadataLookupResponse)
      extends MetadataServiceResponse
  final case class CostFailure(id: WorkflowId, reason: Throwable) extends MetadataServiceFailure

  final case class MetadataWriteSuccess(events: Iterable[MetadataEvent]) extends MetadataServiceResponse
  final case class MetadataWriteFailure(reason: Throwable, events: Iterable[MetadataEvent])
      extends MetadataServiceFailure

  sealed abstract class WorkflowValidationResponse extends MetadataServiceResponse
  case object RecognizedWorkflowId extends WorkflowValidationResponse
  case object UnrecognizedWorkflowId extends WorkflowValidationResponse
  final case class FailedToCheckWorkflowId(cause: Throwable) extends WorkflowValidationResponse

  sealed abstract class FetchWorkflowArchiveStatusAndEndTimeResponse extends MetadataServiceResponse
  final case class WorkflowMetadataArchivedStatusAndEndTime(archiveStatus: MetadataArchiveStatus,
                                                            endTime: Option[OffsetDateTime]
  ) extends FetchWorkflowArchiveStatusAndEndTimeResponse
  final case class FailedToGetArchiveStatusAndEndTime(reason: Throwable)
      extends FetchWorkflowArchiveStatusAndEndTimeResponse

  sealed abstract class MetadataQueryResponse extends MetadataServiceResponse
  final case class WorkflowQuerySuccess(response: WorkflowQueryResponse, meta: Option[QueryMetadata])
      extends MetadataQueryResponse
  final case class WorkflowQueryFailure(reason: Throwable) extends MetadataQueryResponse

  implicit private class EnhancedWomTraversable(val womValues: Iterable[WomValue]) extends AnyVal {
    def toEvents(metadataKey: MetadataKey, fileMap: Map[Path, Path]): List[MetadataEvent] = if (womValues.isEmpty) {
      List(MetadataEvent.empty(metadataKey.copy(key = s"${metadataKey.key}[]")))
    } else {
      womValues.toList.zipWithIndex
        .flatMap { case (value, index) =>
          womValueToMetadataEvents(metadataKey.copy(key = s"${metadataKey.key}[$index]"), value, fileMap)
        }
    }
  }

  /**
    * @param metadataKey Coordinates uniquely identifying what this event describes, such as `outputs:my_workflow.file_array[5]`.
    * @param womValue The value to be serialized into event(s). Complex types are flattened into multiple individual events.
    * @param fileMap Files that copy/move after the workflow completes need to have their new location in metadata.
    *                When creating events for files only, check the map for a possible new location, or empty map for no-op.
    * @return Flat list of metadata events ready to be written to the database.
    */
  def womValueToMetadataEvents(metadataKey: MetadataKey,
                               womValue: WomValue,
                               fileMap: Map[Path, Path] = Map.empty
  ): Iterable[MetadataEvent] = womValue match {
    case WomArray(_, valueSeq) => valueSeq.toEvents(metadataKey, fileMap)
    case WomMap(_, valueMap) =>
      if (valueMap.isEmpty) {
        List(MetadataEvent.empty(metadataKey))
      } else {
        valueMap.toList flatMap { case (key, value) =>
          womValueToMetadataEvents(metadataKey.copy(key = metadataKey.key + s":${key.valueString}"), value, fileMap)
        }
      }
    case objectLike: WomObjectLike =>
      if (objectLike.values.isEmpty) {
        List(MetadataEvent.empty(metadataKey))
      } else {
        objectLike.values.toList flatMap { case (key, value) =>
          womValueToMetadataEvents(metadataKey.copy(key = metadataKey.key + s":$key"), value, fileMap)
        }
      }
    case WomOptionalValue(_, Some(value)) =>
      womValueToMetadataEvents(metadataKey, value, fileMap)
    case WomPair(left, right) =>
      womValueToMetadataEvents(metadataKey.copy(key = metadataKey.key + ":left"), left, fileMap) ++
        womValueToMetadataEvents(metadataKey.copy(key = metadataKey.key + ":right"), right, fileMap)
    case file: WomSingleFile =>
      // Our lookup key is a string; to avoid exceptions, stringify paths instead of pathifying the string
      val stringifiedMap: Map[String, String] = fileMap map { case (src: Path, dst: Path) =>
        src.pathAsString -> dst.pathAsString
      }
      // Why? When we copy/move final outputs, we need to map the original file to the destination file.
      val mappedFile: WomSingleFile = stringifiedMap.get(file.valueString) match {
        case Some(dst) =>
          WomSingleFile(dst)
        case None =>
          file
      }
      List(MetadataEvent(metadataKey, MetadataValue(mappedFile)))
    case value =>
      List(MetadataEvent(metadataKey, MetadataValue(value)))
  }

  def throwableToMetadataEvents(metadataKey: MetadataKey,
                                t: Throwable,
                                failureIndex: Int = Random.nextInt(Int.MaxValue)
  ): List[MetadataEvent] = {
    val emptyCauseList = List(
      MetadataEvent.empty(metadataKey.copy(key = metadataKey.key + s"[$failureIndex]:causedBy[]"))
    )
    val metadataKeyAndFailureIndex = s"${metadataKey.key}[$failureIndex]"

    t match {
      case aggregation: ThrowableAggregation =>
        val message = List(
          MetadataEvent(metadataKey.copy(key = s"$metadataKeyAndFailureIndex:message"),
                        MetadataValue(aggregation.exceptionContext)
          )
        )
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
        val message = List(
          MetadataEvent(metadataKey.copy(key = s"$metadataKeyAndFailureIndex:message"),
                        MetadataValue(aggregation.exceptionContext)
          )
        )
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
        val message = List(
          MetadataEvent(metadataKey.copy(key = s"$metadataKeyAndFailureIndex:message"), MetadataValue(t.getMessage))
        )
        val causeKey = metadataKey.copy(key = s"$metadataKeyAndFailureIndex:causedBy")
        val cause = Option(t.getCause) map { cause =>
          throwableToMetadataEvents(causeKey, cause, 0)
        } getOrElse emptyCauseList
        message ++ cause
    }
  }
}
