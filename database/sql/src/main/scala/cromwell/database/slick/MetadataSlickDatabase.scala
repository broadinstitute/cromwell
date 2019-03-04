package cromwell.database.slick

import java.sql.Timestamp
import java.util.concurrent.atomic.AtomicInteger

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

  def NOT(action: SQLActionBuilder): SQLActionBuilder = concat(sql"NOT ", action)

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

    val summaryTableAlias = "summaryTable"
    val labelsTableCustomName = "labelTable"

    val openingStatement =
      sql"""select #$summaryTableAlias.WORKFLOW_EXECUTION_UUID, #$summaryTableAlias.WORKFLOW_NAME, #$summaryTableAlias.WORKFLOW_STATUS, #$summaryTableAlias.START_TIMESTAMP, #$summaryTableAlias.END_TIMESTAMP, #$summaryTableAlias.SUBMISSION_TIMESTAMP, #$summaryTableAlias.WORKFLOW_METADATA_SUMMARY_ENTRY_ID
            | FROM WORKFLOW_METADATA_SUMMARY_ENTRY #$summaryTableAlias, CUSTOM_LABEL_ENTRY #$labelsTableCustomName
            | WHERE #$summaryTableAlias.workflow_execution_uuid = #$labelsTableCustomName.workflow_execution_uuid""".stripMargin

    val statusConstraint = NonEmptyList.fromList(workflowStatuses.toList.map(status => sql"""#$summaryTableAlias.WORKFLOW_STATUS=$status""")).map(OR).toList
    val nameConstraint = NonEmptyList.fromList(workflowNames.toList.map(name => sql"""#$summaryTableAlias.WORKFLOW_NAME=$name""")).map(OR).toList
    val idConstraint = NonEmptyList.fromList(workflowExecutionUuids.toList.map(uuid => sql"""#$summaryTableAlias.WORKFLOW_EXECUTION_UUID=$uuid""")).map(OR).toList
    val submissionTimeConstraint = submissionTimestampOption.map(ts => sql"""#$summaryTableAlias.SUBMISSION_TIMESTAMP>=$ts""").toList
    val startTimeConstraint = startTimestampOption.map(ts => sql"""#$summaryTableAlias.START_TIMESTAMP>=$ts""").toList
    val endTimeConstraint = endTimestampOption.map(ts => sql"""#$summaryTableAlias.END_TIMESTAMP<=$ts""").toList

    val mixinTableCounter = new AtomicInteger(0)

    // *ALL* of the labelAnd list of KV pairs must exist:
    def labelExists(labelKey: String, labelValue: String) = {
      val tableName = s"labelsMixin" + mixinTableCounter.getAndIncrement()
      sql"""EXISTS(SELECT * from CUSTOM_LABEL_ENTRY #$tableName WHERE ((#$tableName.WORKFLOW_EXECUTION_UUID = #$summaryTableAlias.WORKFLOW_EXECUTION_UUID) AND (#$tableName.CUSTOM_LABEL_KEY = $labelKey) AND (#$tableName.CUSTOM_LABEL_VALUE = $labelValue)))"""
    }
    val labelsAndConstraint = NonEmptyList.fromList(labelAndKeyLabelValues.toList.map { case (labelKey, labelValue) => labelExists(labelKey, labelValue) } ).map(AND).toList

    // At least one of the labelOr list of KV pairs must exist:
    val labelOrConstraint = NonEmptyList.fromList(labelOrKeyLabelValues.toList.map { case (k, v) =>
      AND(NonEmptyList.of(sql"#$labelsTableCustomName.custom_label_key=$k") :+ sql"#$labelsTableCustomName.custom_label_value=$v")
    }).map(OR).toList

    // *ALL* of the excludeLabelOr list of KV pairs must *NOT* exist:
    val excludeLabelsOrConstraint = NonEmptyList.fromList(excludeLabelOrValues.toList.map { case (labelKey, labelValue) => NOT(labelExists(labelKey, labelValue)) } ).map(AND).toList

    // At least one of the excludeLabelAnd list of KV pairs must *NOT* exist:
    val excludeLabelsAndConstraint = NonEmptyList.fromList(excludeLabelAndValues.toList.map { case (labelKey, labelValue) => NOT(labelExists(labelKey, labelValue)) } ).map(OR).toList

    val includeSubworkflowsConstraint = if (includeSubworkflows) List.empty else {
      val tableName = s"subworkflowMixin" + mixinTableCounter.getAndIncrement()
      List(sql"""NOT EXISTS(SELECT * from METADATA_ENTRY #$tableName WHERE ((#$tableName.WORKFLOW_EXECUTION_UUID = #$summaryTableAlias.WORKFLOW_EXECUTION_UUID) AND (#$tableName.METADATA_KEY = 'parentWorkflowId') AND (#$tableName.METADATA_VALUE IS NOT NULL)))""")
    }

    val constraintSet = AND(
      NonEmptyList.of(openingStatement) ++
      statusConstraint ++
      nameConstraint ++
      idConstraint ++
      submissionTimeConstraint ++
      startTimeConstraint ++
      endTimeConstraint ++
      labelOrConstraint ++
      labelsAndConstraint ++
      excludeLabelsOrConstraint ++
      excludeLabelsAndConstraint ++
      includeSubworkflowsConstraint
    )

    val withPaginationCodicil: SQLActionBuilder = (page, pageSize) match {
      case (Some(p), Some(ps)) => concat(constraintSet, sql""" LIMIT ${Integer.max(p-1, 0) * ps},$ps """)
      case (None, Some(ps)) => concat(constraintSet, sql""" LIMIT 0,$ps """)
      case _ => constraintSet
    }

    val result = withPaginationCodicil.as[(String, Option[String], Option[String], Option[Timestamp], Option[Timestamp], Option[Timestamp], Option[Long])]

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
