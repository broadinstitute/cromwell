package cromwell.database.slick

import java.sql.Timestamp

import cromwell.database.sql.MetadataSqlDatabase
import cromwell.database.sql.tables.{Metadatum, WorkflowMetadataSummary}

import scala.concurrent.{ExecutionContext, Future}
import scalaz._

trait MetadataSlickDatabase extends MetadataSqlDatabase {
  this: SlickDatabase with SummarizingSlickDatabase =>

  import dataAccess.driver.api._

  override def addMetadata(events: Iterable[Metadatum])
                          (implicit ec: ExecutionContext): Future[Unit] = {
    val action = dataAccess.metadataAutoInc ++= events
    runTransaction(action).map(_ => ())
  }

  override def queryMetadataEvents(workflowUuid: String)
                                  (implicit ec: ExecutionContext): Future[Seq[Metadatum]] = {
    val action = dataAccess.metadataByWorkflowUuid(workflowUuid).result
    runTransaction(action)
  }

  override def queryMetadataEvents(workflowUuid: String,
                                   key: String)
                                  (implicit ec: ExecutionContext): Future[Seq[Metadatum]] = {
    val action = dataAccess.metadataByWorkflowUuidAndKey(workflowUuid, key).result
    runTransaction(action)
  }

  override def queryMetadataEvents(workflowUuid: String,
                                   callFqn: String,
                                   index: Option[Int],
                                   attempt: Int)
                                  (implicit ec: ExecutionContext): Future[Seq[Metadatum]] = {
    val action = dataAccess.metadataByWorkflowUuidAndCallFqnAndIndexAndAttempt(workflowUuid, callFqn, index,
      attempt).result
    runTransaction(action)
  }

  override def queryMetadataEvents(workflowUuid: String,
                                   key: String,
                                   callFqn: String,
                                   index: Option[Int],
                                   attempt: Int)
                                  (implicit ec: ExecutionContext): Future[Seq[Metadatum]] = {
    val action = dataAccess.metadataByWorkflowUuidAndKeyAndCallFqnAndIndexAndAttempt(workflowUuid, key, callFqn, index,
      attempt).result
    runTransaction(action)
  }

  override def queryMetadataEventsWithWildcardKeys(workflowUuid: String,
                                                   wildcardKeys: NonEmptyList[String],
                                                   requireEmptyJobKey: Boolean)
                                                  (implicit ec: ExecutionContext): Future[Seq[Metadatum]] = {
    val action = dataAccess.queryMetadataMatchingAnyWildcardKeys(workflowUuid, wildcardKeys, requireEmptyJobKey).result
    runTransaction(action)
  }

  override def queryMetadataEventsWithoutWildcardKeys(workflowUuid: String,
                                                      wildcardKeys: NonEmptyList[String],
                                                      requireEmptyJobKey: Boolean)
                                                     (implicit ec: ExecutionContext): Future[Seq[Metadatum]] = {
    val action = dataAccess.queryMetadataNotMatchingAnyWildcardKeys(workflowUuid, wildcardKeys,
      requireEmptyJobKey).result
    runTransaction(action)
  }

  private def updateMetadata(buildUpdatedSummary:
                             (Option[WorkflowMetadataSummary], Seq[Metadatum]) => WorkflowMetadataSummary)
                            (metadataByUuid: (String, Seq[Metadatum]))(implicit ec: ExecutionContext): DBIO[Unit] = {
    val (uuid, metadataForUuid) = metadataByUuid
    for {
    // There might not be a preexisting summary for a given UUID, so `headOption` the result
      existingSummary <- dataAccess.workflowMetadataSummariesByUuid(uuid).result.headOption
      updatedSummary = buildUpdatedSummary(existingSummary, metadataForUuid)
      _ <- upsertWorkflowMetadataSummary(updatedSummary)
    } yield ()
  }

  private def upsertWorkflowMetadataSummary(summary: WorkflowMetadataSummary)
                                           (implicit ec: ExecutionContext): DBIO[Unit] = {
    if (useSlickUpserts) {
      for {
        _ <- dataAccess.workflowMetadataSummaryAutoInc.insertOrUpdate(summary)
      } yield ()
    } else {
      for {
        updateCount <- dataAccess.workflowMetadataSummariesByUuid(summary.workflowUuid).update(summary)
        _ <- updateCount match {
          case 0 => dataAccess.workflowMetadataSummaryAutoInc += summary
          case 1 => DBIO.successful(Unit)
          case _ => DBIO.failed(new RuntimeException(s"Unexpected summary update count $updateCount"))
        }
      } yield ()
    }
  }

  def refreshMetadataSummaries(startTimestamp: Option[Timestamp],
                               key1: String, key2: String, key3: String, key4: String,
                               buildUpdatedSummary: (Option[WorkflowMetadataSummary], Seq[Metadatum]) =>
                                 WorkflowMetadataSummary)
                              (implicit ec: ExecutionContext): Future[Long] = {
    val action = for {
      startIdOption <- getSummaryStatusMaximumId("WORKFLOW_METADATA_SUMMARY", "METADATA_JOURNAL")
      startId = startIdOption.getOrElse(0L) + 1L
      metadata <- dataAccess.metadataWithIdAndTimestampGreaterThanOrEqual(startId, startTimestamp,
        key1, key2, key3, key4).result
      metadataByWorkflowUuid = metadata.groupBy(_.workflowUuid)
      _ <- DBIO.sequence(metadataByWorkflowUuid map updateMetadata(buildUpdatedSummary))
      maximumId = maximumOrZero(metadata.map(_.metadatumId.get))
      _ <- upsertSummaryStatusMaximumId("WORKFLOW_METADATA_SUMMARY", "METADATA_JOURNAL", maximumId)
    } yield maximumId

    runTransaction(action)
  }

  def getStatus(workflowUuid: String)(implicit ec: ExecutionContext): Future[Option[String]] = {
    val action = dataAccess.workflowStatusByUuid(workflowUuid).result.headOption
    // The workflow might not exist, so `headOption`.  But even if the workflow does exist, the status might be None.
    // So flatten the Option[Option[String]] to Option[String].
    runTransaction(action).map(_.flatten)
  }

  def workflowExists(possibleWorkflowId: String)(implicit ec: ExecutionContext): Future[Boolean] = {
    val action = dataAccess.workflowExists(possibleWorkflowId).result
    runTransaction(action)
  }

  override def queryWorkflowSummaries(statuses: Set[String], names: Set[String], uuids: Set[String],
                                      startDate: Option[Timestamp], endDate: Option[Timestamp],
                                      page: Option[Int], pageSize: Option[Int])
                                     (implicit ec: ExecutionContext): Future[Traversable[WorkflowMetadataSummary]] = {
    val action = dataAccess.queryWorkflowSummaries(statuses, names, uuids, startDate, endDate, page, pageSize).result
    runTransaction(action)
  }

  override def countWorkflowSummaries(statuses: Set[String], names: Set[String], uuids: Set[String],
                                      startDate: Option[Timestamp], endDate: Option[Timestamp])
                                     (implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.countWorkflowSummaries(statuses, names, uuids, startDate, endDate).result
    runTransaction(action)
  }
}
