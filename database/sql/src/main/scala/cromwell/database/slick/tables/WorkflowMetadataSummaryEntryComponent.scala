package cromwell.database.slick.tables

import java.sql.Timestamp

import cromwell.database.sql.tables.WorkflowMetadataSummaryEntry

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

  val workflowStatusesForWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowMetadataSummaryEntry <- workflowMetadataSummaryEntries
      if workflowMetadataSummaryEntry.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowMetadataSummaryEntry.workflowStatus
  )

  def filterWorkflowMetadataSummaryEntries(parentIdWorkflowMetadataKey: String,
                                           workflowStatuses: Set[String],
                                           workflowNames: Set[String],
                                           workflowExecutionUuids: Set[String],
                                           labelAndKeyValues: Set[(String, String)],
                                           labelOrKeyValues: Set[(String, String)],
                                           excludeLabelAndValues: Set[(String,String)],
                                           excludeLabelOrValues: Set[(String,String)],
                                           submissionTimestampOption: Option[Timestamp],
                                           startTimestampOption: Option[Timestamp],
                                           endTimestampOption: Option[Timestamp],
                                           includeSubworkflows: Boolean): WorkflowMetadataSummaryEntries => Rep[Boolean] = {
    val include: Rep[Boolean] = true
    val exclude: Rep[Boolean] = false

    { workflowMetadataSummaryEntry =>
      // All query parameters are either Options or Sets, so they might have no values specified at all.  The general
      // pattern for these criteria is to map Options and map/reduceLeftOption Sets, resulting in optional filters.

      // All fields but UUID are nullable, necessitating the folds.  If these fields are null in the database we want to
      // filter the row if the relevant filter has been specified.
      val startTimestampFilter = startTimestampOption.
        map(startTimestamp => workflowMetadataSummaryEntry.startTimestamp.fold(ifEmpty = exclude)(_ >= startTimestamp))
      val endTimestampFilter = endTimestampOption.
        map(endTimestamp => workflowMetadataSummaryEntry.endTimestamp.fold(ifEmpty = exclude)(_ <= endTimestamp))
      val submissionTimestampFilter = submissionTimestampOption.
        map(submissionTimestamp => workflowMetadataSummaryEntry.submissionTimestamp.fold(ifEmpty = exclude)(_ >= submissionTimestamp))
      // Names, UUIDs, and statuses are potentially multi-valued, the reduceLeftOption ORs together any name, UUID, or
      // status criteria to include all matching names, UUIDs, and statuses.
      val workflowNameFilter = workflowNames.
        map(workflowName => workflowMetadataSummaryEntry.workflowName.fold(ifEmpty = exclude)(_ === workflowName)).
        reduceLeftOption(_ || _)
      val workflowExecutionUuidFilter = workflowExecutionUuids.
        map(workflowExecutionUuid => workflowMetadataSummaryEntry.workflowExecutionUuid === workflowExecutionUuid).
        reduceLeftOption(_ || _)
      val workflowStatusFilter = workflowStatuses.
        map(workflowStatus => workflowMetadataSummaryEntry.workflowStatus.fold(ifEmpty = exclude)(_ === workflowStatus)).
        reduceLeftOption(_ || _)
      val labelsAndFilter = existsWorkflowLabels(workflowMetadataSummaryEntry, labelAndKeyValues, _ && _)
      val labelsOrFilter = existsWorkflowLabels(workflowMetadataSummaryEntry, labelOrKeyValues, _ || _)
      val excludeLabelsAndFilter = existsWorkflowLabels(workflowMetadataSummaryEntry, excludeLabelAndValues, _ && _).map(v => !v)
      val excludeLabelsOrFilter = existsWorkflowLabels(workflowMetadataSummaryEntry, excludeLabelOrValues, _ || _).map(v => !v)
      val notASubworkflowFilter: Option[Rep[Boolean]] =
        if (includeSubworkflows) None else Some(!metadataEntryExistsForWorkflowExecutionUuid(workflowMetadataSummaryEntry.workflowExecutionUuid, parentIdWorkflowMetadataKey))

      // Put all the optional filters above together in one place.
      val optionalFilters: List[Option[Rep[Boolean]]] = List(
        workflowNameFilter,
        workflowExecutionUuidFilter,
        workflowStatusFilter,
        labelsAndFilter,
        labelsOrFilter,
        excludeLabelsAndFilter,
        excludeLabelsOrFilter,
        startTimestampFilter,
        endTimestampFilter,
        submissionTimestampFilter,
        notASubworkflowFilter
      )
      // Unwrap the optional filters.  If any of these filters are not defined, replace with `include` to include all
      // rows which might otherwise have been filtered.
      val filters = optionalFilters.map(_.getOrElse(include))
      // AND together the filters.  If there are no filters at all return `include`.
      filters.reduceLeftOption(_ && _).getOrElse(include)
    }
  }

  /**
    * Filters workflows associated with the labels, combining the label matches using using op.
    *
    * The op can OR together or AND together the filtered rows.
    */
  private def existsWorkflowLabels(workflowMetadataSummaryEntry: WorkflowMetadataSummaryEntries,
                                   labelKeyValues: Set[(String, String)],
                                   op: (Rep[Boolean], Rep[Boolean]) => Rep[Boolean]): Option[Rep[Boolean]] = {
    labelKeyValues
      .map({
        case (labelKey, labelValue) =>
          existsWorkflowIdLabelKeyAndValue(workflowMetadataSummaryEntry.workflowExecutionUuid, labelKey, labelValue)
      })
      .reduceLeftOption(op)
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
    val filter = filterWorkflowMetadataSummaryEntries(
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
      includeSubworkflows
    )
    workflowMetadataSummaryEntries.filter(filter).length
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
    val filter = filterWorkflowMetadataSummaryEntries(
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
      includeSubworkflows
    )
    val query = workflowMetadataSummaryEntries.filter(filter)
    (page, pageSize) match {
      case (Some(p), Some(ps)) => query.sortBy(_.workflowMetadataSummaryEntryId.desc).drop((p - 1) * ps).take(ps)
      case (None, Some(ps)) => query.take(ps)
      case _ => query
    }
  }
}
