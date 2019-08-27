package cromwell.database.slick.tables

import java.sql.Timestamp

import cats.data
import cats.data.NonEmptyList
import cromwell.database.sql.tables.WorkflowMetadataSummaryEntry
import shapeless.syntax.std.tuple._
import slick.jdbc.{GetResult, PositionedParameters, SQLActionBuilder}

//noinspection SqlDialectInspection
trait WorkflowMetadataSummaryEntryComponent {

  this: DriverComponent with CustomLabelEntryComponent with MetadataEntryComponent =>

  import driver.api.TupleMethods._
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

    def baseProjection = (workflowExecutionUuid, workflowName, workflowStatus, startTimestamp, endTimestamp,
      submissionTimestamp, parentWorkflowExecutionUuid, rootWorkflowExecutionUuid)

    override def * = baseProjection ~ workflowMetadataSummaryEntryId.? <> (WorkflowMetadataSummaryEntry.tupled, WorkflowMetadataSummaryEntry.unapply)

    def forUpdate = baseProjection.shaped <> (
      tuple => WorkflowMetadataSummaryEntry.tupled(tuple :+ None),
      WorkflowMetadataSummaryEntry.unapply(_: WorkflowMetadataSummaryEntry).map(_.reverse.tail.reverse)
    )

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
    } yield workflowMetadataSummaryEntry.forUpdate)

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

    val customLabelEntryTable = quoted("CUSTOM_LABEL_ENTRY")
    val workflowMetadataSummaryEntryTable = quoted("WORKFLOW_METADATA_SUMMARY_ENTRY")

    val workflowExecutionUuidColumn = quoted("WORKFLOW_EXECUTION_UUID")
    val customLabelKeyColumn = quoted("CUSTOM_LABEL_KEY")
    val customLabelValueColumn = quoted("CUSTOM_LABEL_VALUE")
    val parentWorkflowExecutionUuidColumn = quoted("PARENT_WORKFLOW_EXECUTION_UUID")

    val summaryTableAlias = quoted("summaryTable")
    val labelsOrTableAlias = quoted("labelsOrMixin")
    val labelsAndTableAliases = labelAndKeyLabelValues.zipWithIndex.map {
      case (labelPair, i) => quoted(s"labelAndTable$i") -> labelPair
    }.toMap

    val selectColumns = List(
      "WORKFLOW_EXECUTION_UUID",
      "WORKFLOW_NAME",
      "WORKFLOW_STATUS",
      "START_TIMESTAMP",
      "END_TIMESTAMP",
      "SUBMISSION_TIMESTAMP",
      "PARENT_WORKFLOW_EXECUTION_UUID",
      "ROOT_WORKFLOW_EXECUTION_UUID",
      "WORKFLOW_METADATA_SUMMARY_ENTRY_ID",
    )
      .map(quoted)
      .mkString(s"$summaryTableAlias.", ", ", "")

    val select = selectOrCount match {
      case Select =>
        sql"""|SELECT #$selectColumns
              |""".stripMargin
      case Count =>
        sql"""|SELECT COUNT(1)
              |""".stripMargin
    }

    val labelOrJoin = if (labelOrKeyLabelValues.nonEmpty) {
      Option(
        sql"""|  JOIN #$customLabelEntryTable #$labelsOrTableAlias
              |    ON #$summaryTableAlias.#$workflowExecutionUuidColumn
              |      = #$labelsOrTableAlias.#$workflowExecutionUuidColumn
              |""".stripMargin)
    } else None

    val labelAndJoins = labelsAndTableAliases.toList.map { case (labelAndTableAlias, _) =>
      sql"""|  JOIN #$customLabelEntryTable #$labelAndTableAlias
            |    ON #$summaryTableAlias.#$workflowExecutionUuidColumn
            |      = #$labelAndTableAlias.#$workflowExecutionUuidColumn
            |""".stripMargin
    }

    val from = concatNel(NonEmptyList.of(
      sql"""|FROM #$workflowMetadataSummaryEntryTable #$summaryTableAlias
            |""".stripMargin) ++ labelOrJoin.toList ++ labelAndJoins)

    def makeSetConstraint(column: String, elements: Set[String]) = {
      val list = elements.toList.map(element => sql"""#$summaryTableAlias.#${quoted(column)} = $element""")
      NonEmptyList.fromList(list).map(or).toList
    }

    def makeTimeConstraint(column: String, comparison: String, elementOption: Option[Timestamp]) = {
      elementOption.map(element => sql"""#$summaryTableAlias.#${quoted(column)} #$comparison $element""").toList
    }

    val statusConstraint = makeSetConstraint("WORKFLOW_STATUS", workflowStatuses)
    val nameConstraint = makeSetConstraint("WORKFLOW_NAME", workflowNames)
    val idConstraint = makeSetConstraint("WORKFLOW_EXECUTION_UUID", workflowExecutionUuids)

    val submissionTimeConstraint = makeTimeConstraint("SUBMISSION_TIMESTAMP", ">=", submissionTimestampOption)
    val startTimeConstraint = makeTimeConstraint("START_TIMESTAMP", ">=", startTimestampOption)
    val endTimeConstraint = makeTimeConstraint("END_TIMESTAMP", "<=", endTimestampOption)

    // *ALL* of the labelAnd list of KV pairs must exist:
    val labelsAndConstraint = NonEmptyList.fromList(labelsAndTableAliases.toList.map {
      case (labelsAndTableAlias, (labelKey, labelValue)) =>
        and(NonEmptyList.of(
          sql"""#$labelsAndTableAlias.#$customLabelKeyColumn = $labelKey""",
          sql"""#$labelsAndTableAlias.#$customLabelValueColumn = $labelValue""",
        ))
    }).map(and).toList

    // At least one of the labelOr list of KV pairs must exist:
    val labelOrConstraint = NonEmptyList.fromList(labelOrKeyLabelValues.toList.map {
      case (labelKey, labelValue) =>
        and(NonEmptyList.of(
          sql"""#$labelsOrTableAlias.#$customLabelKeyColumn = $labelKey""",
          sql"""#$labelsOrTableAlias.#$customLabelValueColumn = $labelValue""",
        ))
    }).map(or).toList


    var mixinTableCounter = 0

    def labelExists(labelKey: String, labelValue: String) = {
      val tableName = quoted(s"labelsMixin" + mixinTableCounter)
      mixinTableCounter += 1
      sql"""|EXISTS (
            |  SELECT 1 FROM #$customLabelEntryTable #$tableName
            |    WHERE (
            |      (#$tableName.#$workflowExecutionUuidColumn = #$summaryTableAlias.#$workflowExecutionUuidColumn)
            |        AND (#$tableName.#$customLabelKeyColumn = $labelKey)
            |        AND (#$tableName.#$customLabelValueColumn = $labelValue)
            |    )
            |)
            |""".stripMargin
    }

    // *ALL* of the excludeLabelOr list of KV pairs must *NOT* exist:
    val excludeLabelsOrConstraint = NonEmptyList.fromList(excludeLabelOrValues.toList map {
      case (labelKey, labelValue) => not(labelExists(labelKey, labelValue))
    }).map(and).toList

    // At least one of the excludeLabelAnd list of KV pairs must *NOT* exist:
    val excludeLabelsAndConstraint = NonEmptyList.fromList(excludeLabelAndValues.toList.map {
      case (labelKey, labelValue) => not(labelExists(labelKey, labelValue))
    }).map(or).toList

    val includeSubworkflowsConstraint = if (includeSubworkflows) List.empty else {
      List(sql"""#$summaryTableAlias.#$parentWorkflowExecutionUuidColumn IS NULL""".stripMargin)
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
      case (Some(p), Some(ps)) => List(sql""" LIMIT #$ps OFFSET #${ps * ((p - 1) max 0)}""")
      case (None, Some(ps)) => List(sql""" LIMIT #$ps OFFSET 0""")
      case _ => List.empty
    }

    val orderByAddendum =
      sql"""|  ORDER BY #${quoted("WORKFLOW_METADATA_SUMMARY_ENTRY_ID")} DESC
            |""".stripMargin

    // NB you can preview the prepared statement created here by using, for example: println(result.statements.head)

    val fullQuery = concatNel(NonEmptyList(mainQuery, orderByAddendum :: paginationAddendum))

    fullQuery.as[WorkflowMetadataSummaryEntry](rconv = GetResult { r =>
      WorkflowMetadataSummaryEntry(r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<)
    })
  }
}
