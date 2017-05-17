package cromwell.database.slick

import java.sql.Timestamp

import cats.data.NonEmptyList
import cromwell.database.sql.MetadataSqlDatabase
import cromwell.database.sql.tables.{CustomLabelEntry, MetadataEntry, WorkflowMetadataSummaryEntry}
import cromwell.database.sql.SqlConverters._

import scala.concurrent.{ExecutionContext, Future}

trait MetadataSlickDatabase extends MetadataSqlDatabase {
  this: SlickDatabase with SummaryStatusSlickDatabase =>

  import dataAccess.driver.api._

  override def addMetadataEntries(metadataEntries: Iterable[MetadataEntry])
                                 (implicit ec: ExecutionContext): Future[Unit] = {
    val action = DBIO.seq(metadataEntries.grouped(insertBatchSize).map(dataAccess.metadataEntries ++= _).toSeq:_*)
    runTransaction(action).map(_ => ())
  }

  override def metadataEntryExists(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Boolean] = {
    val action = dataAccess.metadataEntryExistsForWorkflowExecutionUuid(workflowExecutionUuid).result
    runTransaction(action)
  }

  override def queryMetadataEntries(workflowExecutionUuid: String)
                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action = dataAccess.metadataEntriesForWorkflowExecutionUuid(workflowExecutionUuid).result
    runTransaction(action)
  }

  override def queryMetadataEntries(workflowExecutionUuid: String,
                                    metadataKey: String)
                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action =
      dataAccess.metadataEntriesForWorkflowExecutionUuidAndMetadataKey((workflowExecutionUuid, metadataKey)).result
    runTransaction(action)
  }

  override def queryMetadataEntries(workflowExecutionUuid: String,
                                    callFullyQualifiedName: String,
                                    jobIndex: Option[Int],
                                    jobAttempt: Int)
                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action = dataAccess.
      metadataEntriesForJobKey((workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt)).result
    runTransaction(action)
  }

  override def queryMetadataEntries(workflowUuid: String,
                                    metadataKey: String,
                                    callFullyQualifiedName: String,
                                    jobIndex: Option[Int],
                                    jobAttempt: Int)
                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action = dataAccess.metadataEntriesForJobKeyAndMetadataKey((
      workflowUuid, metadataKey, callFullyQualifiedName, jobIndex, jobAttempt)).result
    runTransaction(action)
  }

  override def queryMetadataEntriesLikeMetadataKeys(workflowExecutionUuid: String,
                                                    metadataKeys: NonEmptyList[String],
                                                    requireEmptyJobKey: Boolean)
                                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action =
      dataAccess.metadataEntriesLikeMetadataKeys(workflowExecutionUuid, metadataKeys, requireEmptyJobKey).result
    runTransaction(action)
  }

  override def queryMetadataEntryNotLikeMetadataKeys(workflowExecutionUuid: String,
                                                     metadataKeys: NonEmptyList[String],
                                                     requireEmptyJobKey: Boolean)
                                                    (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action =
      dataAccess.metadataEntriesNotLikeMetadataKeys(workflowExecutionUuid, metadataKeys, requireEmptyJobKey).result
    runTransaction(action)
  }

  private def updateWorkflowMetadataSummaryEntry(buildUpdatedWorkflowMetadataSummaryEntry:
                                                 (Option[WorkflowMetadataSummaryEntry], Seq[MetadataEntry]) =>
                                                   WorkflowMetadataSummaryEntry)
                                                (workflowExecutionUuuidAndMetadataEntries: (String, Seq[MetadataEntry]))
                                                (implicit ec: ExecutionContext): DBIO[Unit] = {
    val (workflowExecutionUuid, metadataEntries) = workflowExecutionUuuidAndMetadataEntries
    for {
    // There might not be a preexisting summary for a given UUID, so `headOption` the result
      existingWorkflowMetadataSummaryEntry <- dataAccess.
        workflowMetadataSummaryEntriesForWorkflowExecutionUuid(workflowExecutionUuid).result.headOption
      updatedWorkflowMetadataSummaryEntry = buildUpdatedWorkflowMetadataSummaryEntry(
        existingWorkflowMetadataSummaryEntry, metadataEntries)
      _ <- upsertWorkflowMetadataSummaryEntry(updatedWorkflowMetadataSummaryEntry)
    } yield ()
  }

  private def upsertCustomLabelEntry(metadataEntry: MetadataEntry)
                                    (implicit ec: ExecutionContext): DBIO[Unit] = {
    //Extracting the label key from the MetadataEntry key
    val labelKey = metadataEntry.metadataKey.split("\\:", 2)(1)
    val labelValue = metadataEntry.metadataValue.toRawString
    val customLabelEntry = CustomLabelEntry(labelKey, labelValue, metadataEntry.workflowExecutionUuid)
    if (useSlickUpserts) {
      for {
        _ <- dataAccess.customLabelEntryIdsAutoInc.insertOrUpdate(customLabelEntry)
      } yield ()
    }
    else {
       for {
         updateCount <- dataAccess.entryByWorkflowUuidLabelKeyLabelValue((
             metadataEntry.workflowExecutionUuid,
           labelKey,
           labelValue)).update(customLabelEntry)
         _ <- updateCount match {
           case 0 => dataAccess.customLabelEntryIdsAutoInc += customLabelEntry
           case _ => assertUpdateCount("upsertCustomLabelEntry", updateCount, 1)
         }
       } yield()
    }
  }

  private def upsertWorkflowMetadataSummaryEntry(workflowMetadataSummaryEntry: WorkflowMetadataSummaryEntry)
                                                (implicit ec: ExecutionContext): DBIO[Unit] = {
    if (useSlickUpserts) {
      for {
        _ <- dataAccess.workflowMetadataSummaryEntryIdsAutoInc.insertOrUpdate(workflowMetadataSummaryEntry)
      } yield ()
    } else {
      for {
        updateCount <- dataAccess.
          workflowMetadataSummaryEntriesForWorkflowExecutionUuid(workflowMetadataSummaryEntry.workflowExecutionUuid).
          update(workflowMetadataSummaryEntry)
        _ <- updateCount match {
          case 0 => dataAccess.workflowMetadataSummaryEntryIdsAutoInc += workflowMetadataSummaryEntry
          case _ => assertUpdateCount("upsertWorkflowMetadataSummaryEntry", updateCount, 1)
        }
      } yield ()
    }
  }

  override def refreshMetadataSummaryEntries(startMetadataKey: String, endMetadataKey: String, nameMetadataKey: String,
                                             statusMetadataKey: String, labelMetadataKey: String, buildUpdatedSummary:
                                             (Option[WorkflowMetadataSummaryEntry], Seq[MetadataEntry]) =>
                                               WorkflowMetadataSummaryEntry)
                                            (implicit ec: ExecutionContext): Future[Long] = {
    val likeMetadataLabelKey = labelMetadataKey + "%"
    val action = for {
      previousMetadataEntryIdOption <- getSummaryStatusEntryMaximumId(
        "WORKFLOW_METADATA_SUMMARY_ENTRY", "METADATA_ENTRY")
      previousMetadataEntryId = previousMetadataEntryIdOption.getOrElse(0L)
      metadataEntries <- dataAccess.metadataEntriesForIdGreaterThanOrEqual((
        previousMetadataEntryId + 1L, startMetadataKey, endMetadataKey, nameMetadataKey, statusMetadataKey, likeMetadataLabelKey)).result
      metadataWithoutLabels = metadataEntries.filterNot(_.metadataKey.contains(labelMetadataKey)).groupBy(_.workflowExecutionUuid)
      customLabelEntries = metadataEntries.filter(_.metadataKey.contains(labelMetadataKey))
      _ <- DBIO.sequence(metadataWithoutLabels map updateWorkflowMetadataSummaryEntry(buildUpdatedSummary))
      _ <- DBIO.sequence(customLabelEntries map upsertCustomLabelEntry)
      maximumMetadataEntryId = previousOrMaximum(previousMetadataEntryId, metadataEntries.map(_.metadataEntryId.get))
      _ <- upsertSummaryStatusEntryMaximumId(
        "WORKFLOW_METADATA_SUMMARY_ENTRY", "METADATA_ENTRY", maximumMetadataEntryId)
    } yield maximumMetadataEntryId

    runTransaction(action)
  }

  override def getWorkflowStatus(workflowExecutionUuid: String)
                                (implicit ec: ExecutionContext): Future[Option[String]] = {
    val action = dataAccess.workflowStatusesForWorkflowExecutionUuid(workflowExecutionUuid).result.headOption
    // The workflow might not exist, so `headOption`.  But even if the workflow does exist, the status might be None.
    // So flatten the Option[Option[String]] to Option[String].
    runTransaction(action).map(_.flatten)
  }

  override def queryWorkflowSummaries(workflowStatuses: Set[String], workflowNames: Set[String],
                                      workflowExecutionUuids: Set[String], labelKeyLabelValues: Set[(String,String)],
                                      startTimestampOption: Option[Timestamp],
                                      endTimestampOption: Option[Timestamp], page: Option[Int], pageSize: Option[Int])
                                     (implicit ec: ExecutionContext): Future[Seq[WorkflowMetadataSummaryEntry]] = {
    val action = dataAccess.queryWorkflowMetadataSummaryEntries(workflowStatuses, workflowNames, workflowExecutionUuids,
      labelKeyLabelValues, startTimestampOption, endTimestampOption, page, pageSize).result
    runTransaction(action)
  }

  override def countWorkflowSummaries(workflowStatuses: Set[String], workflowNames: Set[String],
                                      workflowExecutionUuids: Set[String], labelKeyLabelValues: Set[(String, String)],
                                      startTimestampOption: Option[Timestamp],
                                      endTimestampOption: Option[Timestamp])
                                     (implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.countWorkflowMetadataSummaryEntries(workflowStatuses, workflowNames, workflowExecutionUuids,
      labelKeyLabelValues, startTimestampOption, endTimestampOption).result
    runTransaction(action)
  }
}
