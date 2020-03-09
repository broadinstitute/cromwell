package cromwell.database.slick

import java.sql.Timestamp

import cats.implicits._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.database.slick.tables.MetadataDataAccessComponent
import cromwell.database.sql.MetadataSqlDatabase
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.joins.{CallOrWorkflowQuery, CallQuery, MetadataJobQueryValue, WorkflowQuery}
import cromwell.database.sql.tables.{CustomLabelEntry, MetadataEntry, WorkflowMetadataSummaryEntry}

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._

object MetadataSlickDatabase {
  def fromParentConfig(parentConfig: Config = ConfigFactory.load): MetadataSlickDatabase = {
    val databaseConfig = SlickDatabase.getDatabaseConfig("metadata", parentConfig)
    new MetadataSlickDatabase(databaseConfig)
  }
}

class MetadataSlickDatabase(originalDatabaseConfig: Config)
  extends SlickDatabase(originalDatabaseConfig)
    with MetadataSqlDatabase
    with SummaryStatusSlickDatabase
    with SummaryQueueSlickDatabase {
  override lazy val dataAccess = new MetadataDataAccessComponent(slickConfig.profile)

  import dataAccess.driver.api._

  override def existsMetadataEntries()(implicit ec: ExecutionContext): Future[Boolean] = {
    val action = dataAccess.metadataEntriesExists.result
    runTransaction(action)
  }

  override def addMetadataEntries(metadataEntries: Iterable[MetadataEntry])
                                 (implicit ec: ExecutionContext): Future[Unit] = {

    if (metadataEntries.isEmpty) Future.successful(()) else {

      val batchesToWrite = metadataEntries.grouped(insertBatchSize).toList
      val insertActions = batchesToWrite.map { batch =>
        val insertMetadata = dataAccess.metadataEntryIdsAutoInc ++= batch
        insertMetadata.flatMap(ids => writeSummaryQueueEntries(ids))
      }
      runTransaction(DBIO.sequence(insertActions)).void
    }
  }

  override def metadataEntryExists(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Boolean] = {
    val action = dataAccess.metadataEntryExistsForWorkflowExecutionUuid(workflowExecutionUuid).result
    runTransaction(action)
  }

  override def metadataSummaryEntryExists(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Boolean] = {
    val action = dataAccess.workflowMetadataSummaryEntryExistsForWorkflowExecutionUuid(workflowExecutionUuid).result
    runTransaction(action)
  }

  override def queryMetadataEntries(workflowExecutionUuid: String,
                                    timeout: Duration)
                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action = dataAccess.metadataEntriesForWorkflowExecutionUuid(workflowExecutionUuid).result
    runTransaction(action, timeout = timeout)
  }

  override def queryMetadataEntries(workflowExecutionUuid: String,
                                    metadataKey: String,
                                    timeout: Duration)
                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action =
      dataAccess.metadataEntriesForWorkflowExecutionUuidAndMetadataKey((workflowExecutionUuid, metadataKey)).result
    runTransaction(action, timeout = timeout)
  }

  override def queryMetadataEntries(workflowExecutionUuid: String,
                                    callFullyQualifiedName: String,
                                    jobIndex: Option[Int],
                                    jobAttempt: Option[Int],
                                    timeout: Duration)
                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action = dataAccess.
      metadataEntriesForJobKey((workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt)).result
    runTransaction(action, timeout = timeout)
  }

  override def queryMetadataEntries(workflowUuid: String,
                                    metadataKey: String,
                                    callFullyQualifiedName: String,
                                    jobIndex: Option[Int],
                                    jobAttempt: Option[Int],
                                    timeout: Duration)
                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action = dataAccess.metadataEntriesForJobKeyAndMetadataKey((
      workflowUuid, metadataKey, callFullyQualifiedName, jobIndex, jobAttempt)).result
    runTransaction(action, timeout = timeout)
  }

  override def queryMetadataEntryWithKeyConstraints(workflowExecutionUuid: String,
                                                    metadataKeysToFilterFor: List[String],
                                                    metadataKeysToFilterOut: List[String],
                                                    metadataJobQueryValue: MetadataJobQueryValue,
                                                    timeout: Duration)
                                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action = metadataJobQueryValue match {
      case CallQuery(callFqn, jobIndex, jobAttempt) =>
        dataAccess.metadataEntriesForJobWithKeyConstraints(workflowExecutionUuid, metadataKeysToFilterFor, metadataKeysToFilterOut, callFqn, jobIndex, jobAttempt).result
      case WorkflowQuery =>
        dataAccess.metadataEntriesWithKeyConstraints(workflowExecutionUuid, metadataKeysToFilterFor, metadataKeysToFilterOut, requireEmptyJobKey = true).result
      case CallOrWorkflowQuery =>
        dataAccess.metadataEntriesWithKeyConstraints(workflowExecutionUuid, metadataKeysToFilterFor, metadataKeysToFilterOut, requireEmptyJobKey = false).result
    }
    runTransaction(action, timeout = timeout)
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

  private def toCustomLabelEntry(metadataEntry: MetadataEntry): CustomLabelEntry = {
    //Extracting the label key from the MetadataEntry key
    val labelKey = metadataEntry.metadataKey.split("\\:", 2)(1)
    val labelValue = metadataEntry.metadataValue.toRawString
    val customLabelEntry = CustomLabelEntry(labelKey, labelValue, metadataEntry.workflowExecutionUuid)
    customLabelEntry
  }

  private def upsertCustomLabelEntry(customLabelEntry: CustomLabelEntry)
                                    (implicit ec: ExecutionContext): DBIO[Unit] = {
    if (useSlickUpserts) {
      for {
        _ <- dataAccess.customLabelEntryIdsAutoInc.insertOrUpdate(customLabelEntry)
      } yield ()
    } else {
      for {
        updateCount <- dataAccess.
          customLabelEntriesForWorkflowExecutionUuidAndLabelKey(
            (customLabelEntry.workflowExecutionUuid, customLabelEntry.customLabelKey)
          ).update(customLabelEntry)
        _ <- updateCount match {
          case 0 => dataAccess.customLabelEntryIdsAutoInc += customLabelEntry
          case _ => assertUpdateCount("upsertCustomLabelEntry", updateCount, 1)
        }
      } yield ()
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

  override def summarizeIncreasing(startMetadataKey: String,
                                   endMetadataKey: String,
                                   nameMetadataKey: String,
                                   statusMetadataKey: String,
                                   submissionMetadataKey: String,
                                   parentWorkflowIdKey: String,
                                   rootWorkflowIdKey: String,
                                   labelMetadataKey: String,
                                   limit: Int,
                                   buildUpdatedSummary:
                                   (Option[WorkflowMetadataSummaryEntry], Seq[MetadataEntry])
                                     => WorkflowMetadataSummaryEntry)
                                  (implicit ec: ExecutionContext): Future[(Long, Long)] = {
    val action = for {
      rawMetadataEntries <- dataAccess.metadataEntriesToSummarizeQuery(limit.toLong).result
      _ <-
        buildMetadataSummaryFromRawMetadataAndWriteToDb(
          rawMetadataEntries = rawMetadataEntries,
          startMetadataKey = startMetadataKey,
          endMetadataKey = endMetadataKey,
          nameMetadataKey = nameMetadataKey,
          statusMetadataKey = statusMetadataKey,
          submissionMetadataKey = submissionMetadataKey,
          parentWorkflowIdKey = parentWorkflowIdKey,
          rootWorkflowIdKey = rootWorkflowIdKey,
          labelMetadataKey = labelMetadataKey,
          buildUpdatedSummary = buildUpdatedSummary
        )
      summarizedMetadataEntryIds = rawMetadataEntries.flatMap(_.metadataEntryId)
      _ <- deleteSummaryQueueEntriesByMetadataJournalIds(summarizedMetadataEntryIds)
      unsummarizedTotal <- countSummaryQueueEntries()
    } yield (summarizedMetadataEntryIds.length.toLong, unsummarizedTotal.toLong)

    runTransaction(action)
  }

  override def summarizeDecreasing(summaryNameDecreasing: String,
                                   summaryNameIncreasing: String,
                                   startMetadataKey: String,
                                   endMetadataKey: String,
                                   nameMetadataKey: String,
                                   statusMetadataKey: String,
                                   submissionMetadataKey: String,
                                   parentWorkflowIdKey: String,
                                   rootWorkflowIdKey: String,
                                   labelMetadataKey: String,
                                   limit: Int,
                                   buildUpdatedSummary:
                                   (Option[WorkflowMetadataSummaryEntry], Seq[MetadataEntry])
                                     => WorkflowMetadataSummaryEntry)
                                  (implicit ec: ExecutionContext): Future[(Long, Long)] = {
    val action = for {
      previousExistingMetadataEntryIdOption <- getSummaryStatusEntrySummaryPosition(summaryNameDecreasing)
      previousInitializedMetadataEntryIdOption <- previousExistingMetadataEntryIdOption match {
        case Some(value) => DBIO.successful(Option(value).filter(_ > 0))
        case None => getSummaryStatusEntrySummaryPosition(summaryNameIncreasing).map(_.map(_ + 1L))
      }
      newMinimumMetadataEntryId <- previousInitializedMetadataEntryIdOption match {
        case None => DBIO.successful(0L)
        case Some(startingValue) =>
          val summaryPosition = 0L max (startingValue - limit)
          for {
            rawMetadataEntries <- dataAccess.metadataEntriesForIdRange((summaryPosition, startingValue - 1L)).result
            _ <-
              buildMetadataSummaryFromRawMetadataAndWriteToDb(
                rawMetadataEntries = rawMetadataEntries,
                startMetadataKey = startMetadataKey,
                endMetadataKey = endMetadataKey,
                nameMetadataKey = nameMetadataKey,
                statusMetadataKey = statusMetadataKey,
                submissionMetadataKey = submissionMetadataKey,
                parentWorkflowIdKey = parentWorkflowIdKey,
                rootWorkflowIdKey = rootWorkflowIdKey,
                labelMetadataKey = labelMetadataKey,
                buildUpdatedSummary = buildUpdatedSummary
              )
            _ <- upsertSummaryStatusEntrySummaryPosition(summaryNameDecreasing, summaryPosition)
          } yield summaryPosition
      }
      rowsProcessed = previousExistingMetadataEntryIdOption.map(_ - newMinimumMetadataEntryId).getOrElse(0L)
    } yield (rowsProcessed, newMinimumMetadataEntryId)

    runTransaction(action)
  }

  private def buildMetadataSummaryFromRawMetadataAndWriteToDb(rawMetadataEntries: Seq[MetadataEntry],
                                                              startMetadataKey: String,
                                                              endMetadataKey: String,
                                                              nameMetadataKey: String,
                                                              statusMetadataKey: String,
                                                              submissionMetadataKey: String,
                                                              parentWorkflowIdKey: String,
                                                              rootWorkflowIdKey: String,
                                                              labelMetadataKey: String,
                                                              buildUpdatedSummary:
                                                              (Option[WorkflowMetadataSummaryEntry], Seq[MetadataEntry]) => WorkflowMetadataSummaryEntry
                                                             )(implicit ec: ExecutionContext): DBIO[Unit] = {

    val exactMatchMetadataKeys = Set(startMetadataKey, endMetadataKey, nameMetadataKey, statusMetadataKey, submissionMetadataKey, parentWorkflowIdKey, rootWorkflowIdKey)
    val startsWithMetadataKeys = Set(labelMetadataKey)

    val metadataEntries = rawMetadataEntries filter { entry =>
      entry.callFullyQualifiedName.isEmpty && entry.jobIndex.isEmpty && entry.jobAttempt.isEmpty &&
        (exactMatchMetadataKeys.contains(entry.metadataKey) || startsWithMetadataKeys.exists(entry.metadataKey.startsWith))
    }
    val metadataWithoutLabels = metadataEntries
      .filterNot(_.metadataKey.contains(labelMetadataKey)) // Why are these "contains" while the filtering is "starts with"?
      .groupBy(_.workflowExecutionUuid)
    val customLabelEntries = metadataEntries.filter(_.metadataKey.contains(labelMetadataKey))

    for {
      _ <- DBIO.sequence(metadataWithoutLabels map updateWorkflowMetadataSummaryEntry(buildUpdatedSummary))
      _ <- DBIO.sequence(customLabelEntries map toCustomLabelEntry map upsertCustomLabelEntry)
    } yield ()
  }

  override def updateMetadataArchiveStatus(workflowExecutionUuid: String, newArchiveStatus: Option[String]): Future[Int] = {
    val action = dataAccess.metadataArchiveStatusByWorkflowIdOrRootWorkflowId(workflowExecutionUuid).update(newArchiveStatus)
    runTransaction(action)
  }

  override def getWorkflowStatus(workflowExecutionUuid: String)
                                (implicit ec: ExecutionContext): Future[Option[String]] = {
    val action = dataAccess.workflowStatusesForWorkflowExecutionUuid(workflowExecutionUuid).result.headOption
    // The workflow might not exist, so `headOption`.  But even if the workflow does exist, the status might be None.
    // So flatten the Option[Option[String]] to Option[String].
    runTransaction(action).map(_.flatten)
  }

  override def getWorkflowLabels(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Map[String, String]] = {
    val action = dataAccess.labelsForWorkflowExecutionUuid(workflowExecutionUuid).result
    runTransaction(action).map(_.toMap)
  }

  override def getRootAndSubworkflowLabels(rootWorkflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Map[String, Map[String, String]]] = {
    val action = dataAccess.labelsForWorkflowAndSubworkflows(rootWorkflowExecutionUuid).result
    // An empty Map of String workflow IDs to an inner Map of label keys to label values.
    // The outer Map has a default value so any request for a workflow ID not already present
    // will return an empty inner Map.
    val zero: Map[String, Map[String, String]] = Map.empty.withDefaultValue(Map.empty)

    runTransaction(action) map { seq =>
      seq.foldLeft(zero) { case (labels, (id, labelKey, labelValue)) =>
        val labelsForId = labels(id)
        labels + (id -> (labelsForId + (labelKey -> labelValue)))
      }
    }
  }

  override def queryWorkflowSummaries(parentIdWorkflowMetadataKey: String,
                                      workflowStatuses: Set[String],
                                      workflowNames: Set[String],
                                      workflowExecutionUuids: Set[String],
                                      labelAndKeyLabelValues: Set[(String,String)],
                                      labelOrKeyLabelValues: Set[(String,String)],
                                      excludeLabelAndValues: Set[(String,String)],
                                      excludeLabelOrValues: Set[(String,String)],
                                      submissionTimestampOption: Option[Timestamp],
                                      startTimestampOption: Option[Timestamp],
                                      endTimestampOption: Option[Timestamp],
                                      metadataArchiveStatus: Set[Option[String]],
                                      includeSubworkflows: Boolean,
                                      minimumSummaryEntryId: Option[Long],
                                      page: Option[Int],
                                      pageSize: Option[Int])
                                     (implicit ec: ExecutionContext): Future[Seq[WorkflowMetadataSummaryEntry]] = {

    val action = dataAccess.queryWorkflowMetadataSummaryEntries(parentIdWorkflowMetadataKey, workflowStatuses, workflowNames, workflowExecutionUuids,
      labelAndKeyLabelValues, labelOrKeyLabelValues, excludeLabelAndValues, excludeLabelOrValues, submissionTimestampOption, startTimestampOption,
      endTimestampOption, metadataArchiveStatus, includeSubworkflows, minimumSummaryEntryId, page, pageSize)
    runTransaction(action)
  }

  override def countWorkflowSummaries(parentIdWorkflowMetadataKey: String,
                                      workflowStatuses: Set[String],
                                      workflowNames: Set[String],
                                      workflowExecutionUuids: Set[String],
                                      labelAndKeyLabelValues: Set[(String,String)],
                                      labelOrKeyLabelValues: Set[(String,String)],
                                      excludeLabelAndValues: Set[(String,String)],
                                      excludeLabelOrValues: Set[(String,String)],
                                      submissionTimestampOption: Option[Timestamp],
                                      startTimestampOption: Option[Timestamp],
                                      endTimestampOption: Option[Timestamp],
                                      metadataArchiveStatus: Set[Option[String]],
                                      includeSubworkflows: Boolean,
                                      minimumSummaryEntryId: Option[Long])
                                     (implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.countWorkflowMetadataSummaryEntries(parentIdWorkflowMetadataKey, workflowStatuses, workflowNames, workflowExecutionUuids,
      labelAndKeyLabelValues, labelOrKeyLabelValues, excludeLabelAndValues, excludeLabelOrValues, submissionTimestampOption, startTimestampOption,
      endTimestampOption, metadataArchiveStatus, includeSubworkflows, minimumSummaryEntryId)
    runTransaction(action)
  }

  override def deleteNonLabelMetadataForWorkflowAndUpdateArchiveStatus(rootWorkflowId: String, newArchiveStatus: Option[String])(implicit ec: ExecutionContext): Future[Int] = {
    runTransaction {
      for {
        numDeleted <- dataAccess.metadataEntriesWithoutLabelsForRootWorkflowId(rootWorkflowId).delete
        _ <- dataAccess.metadataArchiveStatusByWorkflowIdOrRootWorkflowId(rootWorkflowId).update(newArchiveStatus)
      } yield numDeleted
    }
  }

  override def isRootWorkflow(rootWorkflowId: String)(implicit ec: ExecutionContext): Future[Option[Boolean]] = {
    runTransaction(
      dataAccess.isRootWorkflow(rootWorkflowId).result.headOption
    )
  }

  override def getRootWorkflowId(workflowId: String)(implicit ec: ExecutionContext): Future[Option[String]] = {
    runAction(
      dataAccess.rootWorkflowId(workflowId).result.headOption
    )
  }

  override def queryRootWorkflowIdsByArchiveStatusAndEndedOnOrBeforeThresholdTimestamp(archiveStatus: Option[String], thresholdTimestamp: Timestamp, batchSize: Long)(implicit ec: ExecutionContext): Future[Seq[String]] = {
    runAction(
      dataAccess.rootWorkflowIdsByArchiveStatusAndEndedOnOrBeforeThresholdTimestamp((archiveStatus, thresholdTimestamp, batchSize)).result
    )
  }

  override def countRootWorkflowIdsByArchiveStatusAndEndedOnOrBeforeThresholdTimestamp(archiveStatus: Option[String], thresholdTimestamp: Timestamp)(implicit ec: ExecutionContext): Future[Int] = {
    runAction(
      dataAccess.countRootWorkflowIdsByArchiveStatusAndEndedOnOrBeforeThresholdTimestamp((archiveStatus, thresholdTimestamp)).result
    )
  }
}
