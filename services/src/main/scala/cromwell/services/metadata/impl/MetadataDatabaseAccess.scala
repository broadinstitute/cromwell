package cromwell.services.metadata.impl

import java.time.OffsetDateTime

import cats.Semigroup
import cats.data.NonEmptyList
import cats.syntax.semigroup._
import cromwell.core.{WorkflowId, WorkflowMetadataKeys, WorkflowState}
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.tables.{MetadataEntry, WorkflowMetadataSummaryEntry}
import cromwell.services.ServicesStore
import cromwell.services.metadata.MetadataService.{QueryMetadata, WorkflowQueryResponse}
import cromwell.services.metadata._

import scala.concurrent.{ExecutionContext, Future}

object MetadataDatabaseAccess {

  private lazy val WorkflowMetadataSummarySemigroup = new Semigroup[WorkflowMetadataSummaryEntry] {
    override def combine(summary1: WorkflowMetadataSummaryEntry,
                         summary2: WorkflowMetadataSummaryEntry): WorkflowMetadataSummaryEntry = {
      // Resolve the status if both `this` and `that` have defined statuses.  This will evaluate to `None`
      // if one or both of the statuses is not defined.
      val resolvedStatus = for {
        thisStatus <- summary1.workflowStatus map WorkflowState.fromString
        thatStatus <- summary2.workflowStatus map WorkflowState.fromString
      } yield (thisStatus |+| thatStatus).toString

      WorkflowMetadataSummaryEntry(
        workflowExecutionUuid = summary1.workflowExecutionUuid,
        workflowName = summary1.workflowName orElse summary2.workflowName,
        workflowStatus = resolvedStatus orElse summary1.workflowStatus orElse summary2.workflowStatus,
        startTimestamp = summary1.startTimestamp orElse summary2.startTimestamp,
        endTimestamp = summary1.endTimestamp orElse summary2.endTimestamp
      )
    }
  }

  def baseSummary(workflowUuid: String) = WorkflowMetadataSummaryEntry(workflowUuid, None, None, None, None, None)

  // If visibility is made `private`, there's a bogus warning about this being unused.
  implicit class MetadatumEnhancer(val metadatum: MetadataEntry) extends AnyVal {
    def toSummary: WorkflowMetadataSummaryEntry = {
      val base = baseSummary(metadatum.workflowExecutionUuid)
      metadatum.metadataKey match {
        case WorkflowMetadataKeys.Name => base.copy(workflowName = metadatum.metadataValue)
        case WorkflowMetadataKeys.Status => base.copy(workflowStatus = metadatum.metadataValue)
        case WorkflowMetadataKeys.StartTime =>
          base.copy(startTimestamp = metadatum.metadataValue map OffsetDateTime.parse map { _.toSystemTimestamp })
        case WorkflowMetadataKeys.EndTime =>
          base.copy(endTimestamp = metadatum.metadataValue map OffsetDateTime.parse map { _.toSystemTimestamp })
      }
    }
  }

  private def buildUpdatedSummary(existingSummary: Option[WorkflowMetadataSummaryEntry],
                                  metadataForUuid: Seq[MetadataEntry]): WorkflowMetadataSummaryEntry = {
    implicit val wmss = WorkflowMetadataSummarySemigroup

    val base = existingSummary.getOrElse(baseSummary(metadataForUuid.head.workflowExecutionUuid))
    metadataForUuid.foldLeft(base) {
      case (metadataSummary, metadatum) => metadataSummary |+| metadatum.toSummary
    }
  }
}

trait MetadataDatabaseAccess {
  this: ServicesStore =>

  // Only for tests
  private[impl] def closeDatabaseInterface() = databaseInterface.close()

  def addMetadataEvents(metadataEvents: Iterable[MetadataEvent])(implicit ec: ExecutionContext): Future[Unit] = {
    val metadata = metadataEvents map { metadataEvent =>
      val key = metadataEvent.key
      val workflowUuid = key.workflowId.id.toString
      val timestamp = metadataEvent.offsetDateTime.toSystemTimestamp
      val value = metadataEvent.value map { _.value }
      val valueType = metadataEvent.value map { _.valueType.typeName }
      val jobKey = key.jobKey map { jk => (jk.callFqn, jk.index, jk.attempt) }
      MetadataEntry(
        workflowUuid, jobKey.map(_._1), jobKey.flatMap(_._2), jobKey.map(_._3), key.key, value, valueType, timestamp)
    }

    databaseInterface.addMetadataEntries(metadata)
  }

  private def metadataToMetadataEvents(workflowId: WorkflowId)(metadata: Seq[MetadataEntry]): Seq[MetadataEvent] = {
    metadata map { m =>
      // If callFullyQualifiedName is non-null then attempt will also be non-null and there is a MetadataJobKey.
      val metadataJobKey: Option[MetadataJobKey] = for {
        callFqn <- m.callFullyQualifiedName
        attempt <- m.jobAttempt
      } yield MetadataJobKey(callFqn, m.jobIndex, attempt)

      val key = MetadataKey(workflowId, metadataJobKey, m.metadataKey)
      val value =  for {
        mValue <- m.metadataValue
        mType <- m.metadataValueType
      } yield MetadataValue(mValue, MetadataType.fromString(mType))

      MetadataEvent(key, value, m.metadataTimestamp.toSystemOffsetDateTime)
    }
  }

  def queryMetadataEvents(query: MetadataQuery)(implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    val uuid = query.workflowId.id.toString

    val futureMetadata: Future[Seq[MetadataEntry]] = query match {
      case MetadataQuery(_, None, None, None, None, _) => databaseInterface.queryMetadataEntries(uuid)
      case MetadataQuery(_, None, Some(key), None, None, _) => databaseInterface.queryMetadataEntries(uuid, key)
      case MetadataQuery(_, Some(jobKey), None, None, None, _) =>
        databaseInterface.queryMetadataEntries(uuid, jobKey.callFqn, jobKey.index, jobKey.attempt)
      case MetadataQuery(_, Some(jobKey), Some(key), None, None, _) =>
        databaseInterface.queryMetadataEntries(uuid, key, jobKey.callFqn, jobKey.index, jobKey.attempt)
      case MetadataQuery(_, None, None, Some(includeKeys), None, _) =>
        databaseInterface.
          queryMetadataEntriesLikeMetadataKeys(uuid, includeKeys.map(_ + "%"), requireEmptyJobKey = false)
      case MetadataQuery(_, None, None, None, Some(excludeKeys), _) =>
        databaseInterface.
          queryMetadataEntryNotLikeMetadataKeys(uuid, excludeKeys.map(_ + "%"), requireEmptyJobKey = false)
      case MetadataQuery(_, None, None, Some(includeKeys), Some(excludeKeys), _) => Future.failed(
        new IllegalArgumentException(
          s"Include/Exclude keys may not be mixed: include = $includeKeys, exclude = $excludeKeys"))
      case invalidQuery => Future.failed(new IllegalArgumentException(
        s"Include/Exclude keys are only supported when querying the workflow, not when querying calls: $invalidQuery"))
    }

    futureMetadata map metadataToMetadataEvents(query.workflowId)
  }

  def queryWorkflowOutputs(id: WorkflowId)
                          (implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    val uuid = id.id.toString
    databaseInterface.queryMetadataEntriesLikeMetadataKeys(
      uuid, NonEmptyList.of(s"${WorkflowMetadataKeys.Outputs}:%"), requireEmptyJobKey = true).
      map(metadataToMetadataEvents(id))
  }

  def queryLogs(id: WorkflowId)
               (implicit ec: ExecutionContext): Future[Seq[MetadataEvent]] = {
    import cromwell.services.metadata.CallMetadataKeys._

    val keys = NonEmptyList.of(Stdout, Stderr, BackendLogsPrefix + ":%")
    databaseInterface.queryMetadataEntriesLikeMetadataKeys(id.id.toString, keys, requireEmptyJobKey = false) map
      metadataToMetadataEvents(id)
  }

  def refreshWorkflowMetadataSummaries()(implicit ec: ExecutionContext): Future[Long] = {
    databaseInterface.refreshMetadataSummaryEntries(WorkflowMetadataKeys.StartTime, WorkflowMetadataKeys.EndTime, WorkflowMetadataKeys.Name,
      WorkflowMetadataKeys.Status, MetadataDatabaseAccess.buildUpdatedSummary)
  }

  def getWorkflowStatus(id: WorkflowId)
                       (implicit ec: ExecutionContext): Future[Option[WorkflowState]] = {
    databaseInterface.getWorkflowStatus(id.toString) map { _ map WorkflowState.fromString }
  }

  def workflowExistsWithId(possibleWorkflowId: String)(implicit ec: ExecutionContext): Future[Boolean] = {
    databaseInterface.metadataEntryExists(possibleWorkflowId)
  }

  def queryWorkflowSummaries(queryParameters: WorkflowQueryParameters)
                            (implicit ec: ExecutionContext): Future[(WorkflowQueryResponse, Option[QueryMetadata])] = {
    val workflowSummaries = databaseInterface.queryWorkflowSummaries(
      queryParameters.statuses, queryParameters.names, queryParameters.ids.map(_.toString),
      queryParameters.startDate.map(_.toSystemTimestamp), queryParameters.endDate.map(_.toSystemTimestamp),
      queryParameters.page, queryParameters.pageSize)

    val workflowSummaryCount = databaseInterface.countWorkflowSummaries(
      queryParameters.statuses, queryParameters.names, queryParameters.ids.map(_.toString),
      queryParameters.startDate.map(_.toSystemTimestamp), queryParameters.endDate.map(_.toSystemTimestamp))

    workflowSummaryCount flatMap { count =>
      workflowSummaries map { workflows =>
        (WorkflowQueryResponse(workflows.toSeq map { workflow =>
          MetadataService.WorkflowQueryResult(
            id = workflow.workflowExecutionUuid,
            name = workflow.workflowName,
            status = workflow.workflowStatus,
            start = workflow.startTimestamp map { _.toSystemOffsetDateTime },
            end = workflow.endTimestamp map { _.toSystemOffsetDateTime })
        }),
          //only return metadata if page is defined
          queryParameters.page map { _ => QueryMetadata(queryParameters.page, queryParameters.pageSize, Option(count)) })
      }
    }
  }
}
