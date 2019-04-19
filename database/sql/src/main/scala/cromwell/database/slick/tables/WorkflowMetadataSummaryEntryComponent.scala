package cromwell.database.slick.tables

import java.sql.Timestamp
import java.util.concurrent.atomic.AtomicInteger

import cats.data
import cats.data.NonEmptyList
import cromwell.database.sql.tables.WorkflowMetadataSummaryEntry
import slick.jdbc.{GetResult, PositionedParameters, SQLActionBuilder}

//noinspection SqlDialectInspection
trait WorkflowMetadataSummaryEntryComponent {

  this: DriverComponent with CustomLabelEntryComponent with MetadataEntryComponent =>

  import driver.api._

  class WorkflowMetadataSummaryEntries(tag: Tag)
    extends Table[WorkflowMetadataSummaryEntry](tag, "WORKFLOW_METADATA_SUMMARY_ENTRY") {
    def workflowMetadataSummaryEntryId = column[Long]("WORKFLOW_METADATA_SUMMARY_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID", O.Length(100))

    def workflowName = column[Option[String]]("WORKFLOW_NAME", O.Length(100))

    def workflowStatus = column[Option[String]]("WORKFLOW_STATUS", O.Length(50))

    def startTimestamp = column[Option[Timestamp]]("START_TIMESTAMP")

    def endTimestamp = column[Option[Timestamp]]("END_TIMESTAMP")

    def submissionTimestamp = column[Option[Timestamp]]("SUBMISSION_TIMESTAMP")

    def parentWorkflowExecutionUuid = column[Option[String]]("PARENT_WORKFLOW_EXECUTION_UUID", O.Length(100))

    def rootWorkflowExecutionUuid = column[Option[String]]("ROOT_WORKFLOW_EXECUTION_UUID", O.Length(100))

    override def * = (workflowExecutionUuid, workflowName, workflowStatus, startTimestamp, endTimestamp,
      submissionTimestamp, parentWorkflowExecutionUuid, rootWorkflowExecutionUuid,
      workflowMetadataSummaryEntryId.?) <> (WorkflowMetadataSummaryEntry.tupled, WorkflowMetadataSummaryEntry.unapply)

    def ucWorkflowMetadataSummaryEntryWeu =
      index("UC_WORKFLOW_METADATA_SUMMARY_ENTRY_WEU", workflowExecutionUuid, unique = true)

    def ixWorkflowMetadataSummaryEntryWn = index("IX_WORKFLOW_METADATA_SUMMARY_ENTRY_WN", workflowName, unique = false)

    def ixWorkflowMetadataSummaryEntryWs =
      index("IX_WORKFLOW_METADATA_SUMMARY_ENTRY_WS", workflowStatus, unique = false)

    def ixWorkflowMetadataSummaryEntryPweu =
      index("IX_WORKFLOW_METADATA_SUMMARY_ENTRY_PWEU", parentWorkflowExecutionUuid, unique = false)

    def ixWorkflowMetadataSummaryEntryRweu =
      index("IX_WORKFLOW_METADATA_SUMMARY_ENTRY_RWEU", rootWorkflowExecutionUuid, unique = false)
  }

  val workflowMetadataSummaryEntries = TableQuery[WorkflowMetadataSummaryEntries]

  val workflowMetadataSummaryEntryIdsAutoInc = workflowMetadataSummaryEntries returning
    workflowMetadataSummaryEntries.map(_.workflowMetadataSummaryEntryId)

  val workflowMetadataSummaryEntriesForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowMetadataSummaryEntry <- workflowMetadataSummaryEntries
      if workflowMetadataSummaryEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowMetadataSummaryEntry)

  val workflowMetadataSummaryEntryExistsForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => (for {
      summaryEntry <- workflowMetadataSummaryEntries
      if summaryEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield summaryEntry).exists
  )

  val workflowStatusesForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowMetadataSummaryEntry <- workflowMetadataSummaryEntries
      if workflowMetadataSummaryEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowMetadataSummaryEntry.workflowStatus
  )

  def concat(a: SQLActionBuilder, b: SQLActionBuilder): SQLActionBuilder = {
    SQLActionBuilder(a.queryParts ++ b.queryParts, (p: Unit, pp: PositionedParameters) => {
        a.unitPConv.apply(p, pp)
        b.unitPConv.apply(p, pp)
    })
  }

  def concatNel(nel: NonEmptyList[SQLActionBuilder]): SQLActionBuilder = nel.tail.foldLeft(nel.head) { (acc, next) => concat(acc, next) }

  def and(list: NonEmptyList[SQLActionBuilder]): SQLActionBuilder = if (list.size == 1) list.head else {
    val fullList = data.NonEmptyList.of(sql"(") ++ list.init.flatMap(x => List(x, sql" AND ")) :+ list.last :+ sql")"
    concatNel(fullList)
  }

  def or(list: NonEmptyList[SQLActionBuilder]): SQLActionBuilder = if (list.size == 1) list.head else {
    val fullList = data.NonEmptyList.of(sql"(") ++ list.init.flatMap(x => List(x, sql" OR ")) :+ list.last :+ sql")"
    concatNel(fullList)
  }

  def not(action: SQLActionBuilder): SQLActionBuilder = concat(sql"NOT ", action)

  sealed trait SelectOrCount
  case object Select extends SelectOrCount
  case object Count extends SelectOrCount

  def buildQueryAction(selectOrCount: SelectOrCount,
                       parentIdWorkflowMetadataKey: String,
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
                       includeSubworkflows: Boolean): SQLActionBuilder = {

    val summaryTableAlias = "summaryTable"
    val labelsOrTableAlias = "labelsOrMixin"
    val labelsAndTableAliases = labelAndKeyLabelValues.zipWithIndex.map { case (labelPair, i) => s"labelAndTable$i" -> labelPair }.toMap

    val select = selectOrCount match {
      case Select =>
        sql"""|SELECT #$summaryTableAlias.WORKFLOW_EXECUTION_UUID,
              |  #$summaryTableAlias.WORKFLOW_NAME,
              |  #$summaryTableAlias.WORKFLOW_STATUS,
              |  #$summaryTableAlias.START_TIMESTAMP,
              |  #$summaryTableAlias.END_TIMESTAMP,
              |  #$summaryTableAlias.SUBMISSION_TIMESTAMP,
              |  #$summaryTableAlias.PARENT_WORKFLOW_EXECUTION_UUID,
              |  #$summaryTableAlias.ROOT_WORKFLOW_EXECUTION_UUID,
              |  #$summaryTableAlias.WORKFLOW_METADATA_SUMMARY_ENTRY_ID
              | """.stripMargin
      case Count =>
        sql"""SELECT COUNT(1)
             | """.stripMargin
    }

    val labelOrJoin = if (labelOrKeyLabelValues.nonEmpty) {
      Option(
        sql""" JOIN CUSTOM_LABEL_ENTRY #$labelsOrTableAlias on #$summaryTableAlias.WORKFLOW_EXECUTION_UUID = #$labelsOrTableAlias.WORKFLOW_EXECUTION_UUID
              | """.stripMargin)
    } else None

    val labelAndJoins = labelsAndTableAliases.toList.map { case (labelAndTableAlias, _) =>
        sql""" JOIN CUSTOM_LABEL_ENTRY #$labelAndTableAlias on #$summaryTableAlias.WORKFLOW_EXECUTION_UUID = #$labelAndTableAlias.WORKFLOW_EXECUTION_UUID
             | """.stripMargin
    }

    val from = concatNel(NonEmptyList.of(
      sql"""FROM WORKFLOW_METADATA_SUMMARY_ENTRY #$summaryTableAlias
           | """.stripMargin) ++ labelOrJoin.toList ++ labelAndJoins )

    val statusConstraint = NonEmptyList.fromList(workflowStatuses.toList.map(status => sql"""#$summaryTableAlias.WORKFLOW_STATUS=$status""")).map(or).toList
    val nameConstraint = NonEmptyList.fromList(workflowNames.toList.map(name => sql"""#$summaryTableAlias.WORKFLOW_NAME=$name""")).map(or).toList
    val idConstraint = NonEmptyList.fromList(workflowExecutionUuids.toList.map(uuid => sql"""#$summaryTableAlias.WORKFLOW_EXECUTION_UUID=$uuid""")).map(or).toList
    val submissionTimeConstraint = submissionTimestampOption.map(ts => sql"""#$summaryTableAlias.SUBMISSION_TIMESTAMP>=$ts""").toList
    val startTimeConstraint = startTimestampOption.map(ts => sql"""#$summaryTableAlias.START_TIMESTAMP>=$ts""").toList
    val endTimeConstraint = endTimestampOption.map(ts => sql"""#$summaryTableAlias.END_TIMESTAMP<=$ts""").toList

    // *ALL* of the labelAnd list of KV pairs must exist:
    val labelsAndConstraint = NonEmptyList.fromList(labelsAndTableAliases.toList.map { case (labelsAndTableAlias, (labelKey, labelValue)) =>
      and(NonEmptyList.of(sql"#$labelsAndTableAlias.custom_label_key=$labelKey") :+ sql"#$labelsAndTableAlias.custom_label_value=$labelValue")
    }).map(and).toList

    // At least one of the labelOr list of KV pairs must exist:
    val labelOrConstraint = NonEmptyList.fromList(labelOrKeyLabelValues.toList.map { case (k, v) =>
      and(NonEmptyList.of(sql"#$labelsOrTableAlias.custom_label_key=$k") :+ sql"#$labelsOrTableAlias.custom_label_value=$v")
    }).map(or).toList


    val mixinTableCounter = new AtomicInteger(0)

    def labelExists(labelKey: String, labelValue: String) = {
      val tableName = s"labelsMixin" + mixinTableCounter.getAndIncrement()
      sql"""EXISTS(SELECT 1 from CUSTOM_LABEL_ENTRY #$tableName WHERE ((#$tableName.WORKFLOW_EXECUTION_UUID = #$summaryTableAlias.WORKFLOW_EXECUTION_UUID) AND (#$tableName.CUSTOM_LABEL_KEY = $labelKey) AND (#$tableName.CUSTOM_LABEL_VALUE = $labelValue)))"""
    }

    // *ALL* of the excludeLabelOr list of KV pairs must *NOT* exist:
    val excludeLabelsOrConstraint = NonEmptyList.fromList(excludeLabelOrValues.toList.map { case (labelKey, labelValue) => not(labelExists(labelKey, labelValue)) } ).map(and).toList

    // At least one of the excludeLabelAnd list of KV pairs must *NOT* exist:
    val excludeLabelsAndConstraint = NonEmptyList.fromList(excludeLabelAndValues.toList.map { case (labelKey, labelValue) => not(labelExists(labelKey, labelValue)) } ).map(or).toList

    val includeSubworkflowsConstraint = if (includeSubworkflows) List.empty else {
      List(sql"""#$summaryTableAlias.PARENT_WORKFLOW_EXECUTION_UUID IS NULL""".stripMargin)
    }

    val constraintList =
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

    val where = NonEmptyList.fromList(constraintList) match {
      case Some(constraints) => List(sql"WHERE ", and(constraints))

      // Is this desirable?
      case None => List.empty
    }

    concatNel((NonEmptyList.of(select) :+ from) ++ where)
  }


  def countWorkflowMetadataSummaryEntries(parentIdWorkflowMetadataKey: String,
                                          workflowStatuses: Set[String], workflowNames: Set[String],
                                          workflowExecutionUuids: Set[String],
                                          labelAndKeyLabelValues: Set[(String,String)],
                                          labelOrKeyLabelValues: Set[(String,String)],
                                          excludeLabelAndValues: Set[(String,String)],
                                          excludeLabelOrValues: Set[(String,String)],
                                          submissionTimestampOption: Option[Timestamp],
                                          startTimestampOption: Option[Timestamp],
                                          endTimestampOption: Option[Timestamp],
                                          includeSubworkflows: Boolean) = {
    buildQueryAction(
      selectOrCount = Count,
      parentIdWorkflowMetadataKey,
      workflowStatuses,
      workflowNames,
      workflowExecutionUuids,
      labelAndKeyLabelValues,
      labelOrKeyLabelValues,
      excludeLabelAndValues,
      excludeLabelOrValues,
      submissionTimestampOption,
      startTimestampOption,
      endTimestampOption,
      includeSubworkflows = includeSubworkflows
    ).as[Int].head
  }

  /**
    * Query workflow execution using the filter criteria encapsulated by the `WorkflowExecutionQueryParameters`.
    */
  def queryWorkflowMetadataSummaryEntries(parentIdWorkflowMetadataKey: String,
                                          workflowStatuses: Set[String], workflowNames: Set[String],
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
                                          pageSize: Option[Int]) = {
    val mainQuery = buildQueryAction(
      selectOrCount = Select,
      parentIdWorkflowMetadataKey,
      workflowStatuses,
      workflowNames,
      workflowExecutionUuids,
      labelAndKeyLabelValues,
      labelOrKeyLabelValues,
      excludeLabelAndValues,
      excludeLabelOrValues,
      submissionTimestampOption,
      startTimestampOption,
      endTimestampOption,
      includeSubworkflows = includeSubworkflows
    )

    val paginationAddendum: List[SQLActionBuilder] = (page, pageSize) match {
      case (Some(p), Some(ps)) => List(sql""" LIMIT #${Integer.max(p-1, 0) * ps},#$ps """)
      case (None, Some(ps)) => List(sql""" LIMIT 0,#$ps """)
      case _ => List.empty
    }

    val orderByAddendum = sql""" ORDER BY WORKFLOW_METADATA_SUMMARY_ENTRY_ID DESC
                               | """.stripMargin

    // NB you can preview the prepared statement created here by using, for example: println(result.statements.head)

    concatNel((NonEmptyList.of(mainQuery) :+ orderByAddendum) ++ paginationAddendum)
      .as[WorkflowMetadataSummaryEntry](rconv = GetResult { r =>
      WorkflowMetadataSummaryEntry(r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<)
    })
  }
}
