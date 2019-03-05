package cromwell.database.slick.tables

import java.sql.Timestamp
import java.util.concurrent.atomic.AtomicInteger

import cats.data
import cats.data.NonEmptyList
import cromwell.database.sql.tables.WorkflowMetadataSummaryEntry
import slick.jdbc.{PositionedParameters, SQLActionBuilder, SetParameter}

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

    override def * = (workflowExecutionUuid, workflowName, workflowStatus, startTimestamp, endTimestamp, submissionTimestamp,
      workflowMetadataSummaryEntryId.?) <> (WorkflowMetadataSummaryEntry.tupled, WorkflowMetadataSummaryEntry.unapply)

    def ucWorkflowMetadataSummaryEntryWeu =
      index("UC_WORKFLOW_METADATA_SUMMARY_ENTRY_WEU", workflowExecutionUuid, unique = true)

    def ixWorkflowMetadataSummaryEntryWn = index("IX_WORKFLOW_METADATA_SUMMARY_ENTRY_WN", workflowName, unique = false)

    def ixWorkflowMetadataSummaryEntryWs =
      index("IX_WORKFLOW_METADATA_SUMMARY_ENTRY_WS", workflowStatus, unique = false)
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

  sealed trait SelectOrCount
  case object SELECT extends SelectOrCount
  case object COUNT extends SelectOrCount

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
    val labelsTableCustomName = "labelsOrMixin"

    val select = selectOrCount match {
      case SELECT =>
        sql"""SELECT #$summaryTableAlias.WORKFLOW_EXECUTION_UUID, #$summaryTableAlias.WORKFLOW_NAME, #$summaryTableAlias.WORKFLOW_STATUS, #$summaryTableAlias.START_TIMESTAMP, #$summaryTableAlias.END_TIMESTAMP, #$summaryTableAlias.SUBMISSION_TIMESTAMP, #$summaryTableAlias.WORKFLOW_METADATA_SUMMARY_ENTRY_ID
             | """.stripMargin
      case COUNT =>
        sql"""SELECT COUNT(1)
             | """.stripMargin
    }

    val from = if (labelOrKeyLabelValues.nonEmpty) {
      sql"""FROM WORKFLOW_METADATA_SUMMARY_ENTRY #$summaryTableAlias, CUSTOM_LABEL_ENTRY #$labelsTableCustomName
           | """.stripMargin
    } else {
      sql"""FROM WORKFLOW_METADATA_SUMMARY_ENTRY #$summaryTableAlias
           | """.stripMargin
    }

    val summaryToLabelsJoin = if (labelOrKeyLabelValues.nonEmpty) {
      List(sql"""#$summaryTableAlias.WORKFLOW_EXECUTION_UUID = #$labelsTableCustomName.WORKFLOW_EXECUTION_UUID""")
    } else List.empty

    val statusConstraint = NonEmptyList.fromList(workflowStatuses.toList.map(status => sql"""#$summaryTableAlias.WORKFLOW_STATUS=$status""")).map(OR).toList
    val nameConstraint = NonEmptyList.fromList(workflowNames.toList.map(name => sql"""#$summaryTableAlias.WORKFLOW_NAME=$name""")).map(OR).toList
    val idConstraint = NonEmptyList.fromList(workflowExecutionUuids.toList.map(uuid => sql"""#$summaryTableAlias.WORKFLOW_EXECUTION_UUID=$uuid""")).map(OR).toList
    val submissionTimeConstraint = submissionTimestampOption.map(ts => sql"""#$summaryTableAlias.SUBMISSION_TIMESTAMP>=$ts""").toList
    val startTimeConstraint = startTimestampOption.map(ts => sql"""#$summaryTableAlias.START_TIMESTAMP>=$ts""").toList
    val endTimeConstraint = endTimestampOption.map(ts => sql"""#$summaryTableAlias.END_TIMESTAMP<=$ts""").toList

    val mixinTableCounter = new AtomicInteger(0)

    def labelExists(labelKey: String, labelValue: String) = {
      val tableName = s"labelsMixin" + mixinTableCounter.getAndIncrement()
      sql"""EXISTS(SELECT 1 from CUSTOM_LABEL_ENTRY #$tableName WHERE ((#$tableName.WORKFLOW_EXECUTION_UUID = #$summaryTableAlias.WORKFLOW_EXECUTION_UUID) AND (#$tableName.CUSTOM_LABEL_KEY = $labelKey) AND (#$tableName.CUSTOM_LABEL_VALUE = $labelValue)))"""
    }

    // *ALL* of the labelAnd list of KV pairs must exist:
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
      List(sql"""NOT EXISTS(SELECT 1 from METADATA_ENTRY #$tableName WHERE ((#$tableName.WORKFLOW_EXECUTION_UUID = #$summaryTableAlias.WORKFLOW_EXECUTION_UUID) AND (#$tableName.METADATA_KEY = 'parentWorkflowId') AND (#$tableName.METADATA_VALUE IS NOT NULL)))""")
    }

    val constraintList =
      summaryToLabelsJoin ++
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
      case Some(constraints) => List(sql"WHERE ", AND(constraints))

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
      selectOrCount = COUNT,
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
    ).as[Int]
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

    // NB you can preview the prepared statement created here by using, for example: println(result.statements.head)

    val mainQuery = buildQueryAction(
      selectOrCount = SELECT,
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

    concatNel((NonEmptyList.of(mainQuery) :+ orderByAddendum) ++ paginationAddendum).as[(String, Option[String], Option[String], Option[Timestamp], Option[Timestamp], Option[Timestamp], Option[Long])]
  }
}
