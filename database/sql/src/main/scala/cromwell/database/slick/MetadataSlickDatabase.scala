package cromwell.database.slick

import java.sql.Timestamp

import cats.data
import cats.data.NonEmptyList
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.database.slick.tables.MetadataDataAccessComponent
import cromwell.database.sql.MetadataSqlDatabase
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.joins.{CallOrWorkflowQuery, CallQuery, MetadataJobQueryValue, WorkflowQuery}
import cromwell.database.sql.tables.{CustomLabelEntry, MetadataEntry, WorkflowMetadataSummaryEntry}
import slick.jdbc.{PositionedParameters, SQLActionBuilder, SetParameter}

import scala.concurrent.{ExecutionContext, Future}

//import slick.driver.H2Driver.api._

object MetadataSlickDatabase {
  def fromParentConfig(parentConfig: Config = ConfigFactory.load): MetadataSlickDatabase = {
    val databaseConfig = SlickDatabase.getDatabaseConfig("metadata", parentConfig)
    new MetadataSlickDatabase(databaseConfig)
  }
}

class MetadataSlickDatabase(originalDatabaseConfig: Config)
  extends SlickDatabase(originalDatabaseConfig)
    with MetadataSqlDatabase
    with SummaryStatusSlickDatabase {
  override lazy val dataAccess = new MetadataDataAccessComponent(slickConfig.profile)

  import dataAccess.driver.api._

  override def existsMetadataEntries()(implicit ec: ExecutionContext): Future[Boolean] = {
    val action = dataAccess.metadataEntriesExists.result
    runTransaction(action)
  }

  override def addMetadataEntries(metadataEntries: Iterable[MetadataEntry])
                                 (implicit ec: ExecutionContext): Future[Unit] = {
    val action = DBIO.seq(metadataEntries.grouped(insertBatchSize).map(dataAccess.metadataEntries ++= _).toSeq:_*)
    runAction(action)
  }

  override def metadataEntryExists(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Boolean] = {
    val action = dataAccess.metadataEntryExistsForWorkflowExecutionUuid(workflowExecutionUuid).result
    runTransaction(action)
  }

  override def metadataSummaryEntryExists(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Boolean] = {
    val action = dataAccess.workflowMetadataSummaryEntryExistsForWorkflowExecutionUuid(workflowExecutionUuid).result
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
                                    jobAttempt: Option[Int])
                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action = dataAccess.
      metadataEntriesForJobKey((workflowExecutionUuid, callFullyQualifiedName, jobIndex, jobAttempt)).result
    runTransaction(action)
  }

  override def queryMetadataEntries(workflowUuid: String,
                                    metadataKey: String,
                                    callFullyQualifiedName: String,
                                    jobIndex: Option[Int],
                                    jobAttempt: Option[Int])
                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action = dataAccess.metadataEntriesForJobKeyAndMetadataKey((
      workflowUuid, metadataKey, callFullyQualifiedName, jobIndex, jobAttempt)).result
    runTransaction(action)
  }

  override def queryMetadataEntriesLikeMetadataKeys(workflowExecutionUuid: String,
                                                    metadataKeys: NonEmptyList[String],
                                                    metadataJobQueryValue: MetadataJobQueryValue)
                                                   (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action = metadataJobQueryValue match {
      case CallQuery(callFqn, jobIndex, jobAttempt) =>
        dataAccess.metadataEntriesLikeMetadataKeysWithJob(workflowExecutionUuid, metadataKeys, callFqn, jobIndex, jobAttempt).result
      case WorkflowQuery => dataAccess.metadataEntriesLikeMetadataKeys(workflowExecutionUuid, metadataKeys, requireEmptyJobKey = true).result
      case CallOrWorkflowQuery => dataAccess.metadataEntriesLikeMetadataKeys(workflowExecutionUuid, metadataKeys, requireEmptyJobKey = false).result
    }
      
    runTransaction(action)
  }

  override def queryMetadataEntryNotLikeMetadataKeys(workflowExecutionUuid: String,
                                                     metadataKeys: NonEmptyList[String],
                                                     metadataJobQueryValue: MetadataJobQueryValue)
                                                    (implicit ec: ExecutionContext): Future[Seq[MetadataEntry]] = {
    val action = metadataJobQueryValue match {
      case CallQuery(callFqn, jobIndex, jobAttempt) =>
        dataAccess.metadataEntriesNotLikeMetadataKeysWithJob(workflowExecutionUuid, metadataKeys, callFqn, jobIndex, jobAttempt).result
      case WorkflowQuery => dataAccess.metadataEntriesNotLikeMetadataKeys(workflowExecutionUuid, metadataKeys, requireEmptyJobKey = true).result
      case CallOrWorkflowQuery => dataAccess.metadataEntriesNotLikeMetadataKeys(workflowExecutionUuid, metadataKeys, requireEmptyJobKey = false).result
    }
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

  override def refreshMetadataSummaryEntries(startMetadataKey: String, endMetadataKey: String, nameMetadataKey: String,
                                             statusMetadataKey: String, labelMetadataKey: String, submissionMetadataKey: String,
                                             buildUpdatedSummary: (Option[WorkflowMetadataSummaryEntry], Seq[MetadataEntry]) =>
                                               WorkflowMetadataSummaryEntry)
                                            (implicit ec: ExecutionContext): Future[Long] = {
    val likeMetadataLabelKey = labelMetadataKey + "%"
    val action = for {
      previousMetadataEntryIdOption <- getSummaryStatusEntryMaximumId(
        "WORKFLOW_METADATA_SUMMARY_ENTRY", "METADATA_ENTRY")
      previousMetadataEntryId = previousMetadataEntryIdOption.getOrElse(-1L)
      metadataEntries <- dataAccess.metadataEntriesForIdGreaterThanOrEqual((
        previousMetadataEntryId + 1L, startMetadataKey, endMetadataKey, nameMetadataKey, statusMetadataKey, likeMetadataLabelKey, submissionMetadataKey)).result
      metadataWithoutLabels = metadataEntries.filterNot(_.metadataKey.contains(labelMetadataKey)).groupBy(_.workflowExecutionUuid)
      customLabelEntries = metadataEntries.filter(_.metadataKey.contains(labelMetadataKey))
      _ <- DBIO.sequence(metadataWithoutLabels map updateWorkflowMetadataSummaryEntry(buildUpdatedSummary))
      _ <- DBIO.sequence(customLabelEntries map toCustomLabelEntry map upsertCustomLabelEntry)
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

  override def getWorkflowLabels(workflowExecutionUuid: String)(implicit ec: ExecutionContext): Future[Map[String, String]] = {
    val action = dataAccess.labelsForWorkflowExecutionUuid(workflowExecutionUuid).result
    runTransaction(action).map(_.toMap)
  }

  def concat(a: SQLActionBuilder, b: SQLActionBuilder): SQLActionBuilder = {
    SQLActionBuilder(a.queryParts ++ b.queryParts, new SetParameter[Unit] {
      def apply(p: Unit, pp: PositionedParameters): Unit = {
        a.unitPConv.apply(p, pp)
        b.unitPConv.apply(p, pp)
      }
    })
  }

  def concatNel(nel: NonEmptyList[SQLActionBuilder]): SQLActionBuilder = nel.tail.foldLeft(nel.head) { (acc, next) => concat(acc, next) }

  def AND(list: NonEmptyList[SQLActionBuilder]): SQLActionBuilder = if (list.size == 1) list.head else {
    val fullList = data.NonEmptyList.of(sql"(") ++ list.init.flatMap(x => List(x, sql" AND ")) :+ list.last :+ sql")"
    concatNel(fullList)
  }

  def OR(list: NonEmptyList[SQLActionBuilder]): SQLActionBuilder = if (list.size == 1) list.head else {
    val fullList = data.NonEmptyList.of(sql"(") ++ list.init.flatMap(x => List(x, sql" OR ")) :+ list.last :+ sql")"
    concatNel(fullList)
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
                                      includeSubworkflows: Boolean,
                                      page: Option[Int],
                                      pageSize: Option[Int])
                                     (implicit ec: ExecutionContext): Future[Seq[WorkflowMetadataSummaryEntry]] = {

    val basis =
      sql"""select se.WORKFLOW_EXECUTION_UUID, se.WORKFLOW_NAME, se.WORKFLOW_STATUS, se.START_TIMESTAMP, se.END_TIMESTAMP, se.SUBMISSION_TIMESTAMP, se.WORKFLOW_METADATA_SUMMARY_ENTRY_ID
            | FROM WORKFLOW_METADATA_SUMMARY_ENTRY se, CUSTOM_LABEL_ENTRY le
            | WHERE se.workflow_execution_uuid = le.workflow_execution_uuid""".stripMargin

    val statusConstraint = NonEmptyList.fromList(workflowStatuses.toList.map(status => sql"""se.WORKFLOW_STATUS=$status""")).map(OR).toList
    val nameConstraint = NonEmptyList.fromList(workflowNames.toList.map(name => sql"""se.WORKFLOW_NAME=$name""")).map(OR).toList
    val idConstraint = NonEmptyList.fromList(workflowExecutionUuids.toList.map(uuid => sql"""se.WORKFLOW_EXECUTION_UUID=$uuid""")).map(OR).toList
    val submissionTimeConstraint = submissionTimestampOption.map(ts => sql"""se.SUBMISSION_TIMESTAMP>=$ts""").toList
    val startTimeConstraint = startTimestampOption.map(ts => sql"""se.START_TIMESTAMP>=$ts""").toList
    val endTimeConstraint = endTimestampOption.map(ts => sql"""se.END_TIMESTAMP<=$ts""").toList

    val labelOrConstraint = NonEmptyList.fromList(labelAndKeyLabelValues.toList.map { case (k, v) =>
      AND(NonEmptyList.of(sql"el.custom_label_key=$k") :+ sql"el.custom_label_value=$v")
    }).map(OR).toList

    val includeSubworkflowsConstraint = if (includeSubworkflows) List.empty else List(
      sql"""NOT EXISTS(SELECT * from METADATA_ENTRY WHERE ((WORKFLOW_EXECUTION_UUID = se.WORKFLOW_EXECUTION_UUID) AND (METADATA_KEY = 'parentWorkflowId') AND (METADATA_VALUE IS NOT NULL)))"""
    )

    val sqlNel = NonEmptyList.of(basis) ++ statusConstraint ++ nameConstraint ++ idConstraint ++ submissionTimeConstraint ++ startTimeConstraint ++ endTimeConstraint ++ labelOrConstraint ++ includeSubworkflowsConstraint

    val result = AND(sqlNel).as[(String, Option[String], Option[String], Option[Timestamp], Option[Timestamp], Option[Timestamp], Option[Long])]

    println(result.statements.head)

    runTransaction(result).map(results => results map { case (a,b,c,d,e,f,g) => WorkflowMetadataSummaryEntry(a,b,c,d,e,f,g) })
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
                                      includeSubworkflows: Boolean)
                                     (implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.countWorkflowMetadataSummaryEntries(parentIdWorkflowMetadataKey, workflowStatuses, workflowNames, workflowExecutionUuids,
      labelAndKeyLabelValues, labelOrKeyLabelValues, excludeLabelAndValues, excludeLabelOrValues, submissionTimestampOption, startTimestampOption, endTimestampOption, includeSubworkflows).result
    runTransaction(action)
  }
}
