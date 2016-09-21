package cromwell.database.slick.tables

import java.sql.Timestamp

import cromwell.database.sql.tables.WorkflowMetadataSummaryEntry

trait WorkflowMetadataSummaryEntryComponent {

  this: DriverComponent =>

  import driver.api._

  class WorkflowMetadataSummaryEntries(tag: Tag)
    extends Table[WorkflowMetadataSummaryEntry](tag, "WORKFLOW_METADATA_SUMMARY_ENTRY") {
    def workflowMetadataSummaryEntryId = column[Long]("WORKFLOW_METADATA_SUMMARY_ENTRY_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")

    def workflowName = column[Option[String]]("WORKFLOW_NAME")

    def workflowStatus = column[Option[String]]("WORKFLOW_STATUS")

    def startTimestamp = column[Option[Timestamp]]("START_TIMESTAMP")

    def endTimestamp = column[Option[Timestamp]]("END_TIMESTAMP")

    override def * = (workflowExecutionUuid, workflowName, workflowStatus, startTimestamp, endTimestamp,
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

  def filterWorkflowMetadataSummaryEntries(workflowStatuses: Set[String], workflowNames: Set[String],
                                           workflowExecutionUuids: Set[String], startTimestampOption: Option[Timestamp],
                                           endTimestampOption: Option[Timestamp]):
  (WorkflowMetadataSummaryEntries) => Rep[Boolean] = {
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

      // Put all the optional filters above together in one place.
      val optionalFilters: List[Option[Rep[Boolean]]] = List(
        workflowNameFilter, workflowExecutionUuidFilter, workflowStatusFilter, startTimestampFilter, endTimestampFilter)
      // Unwrap the optional filters.  If any of these filters are not defined, replace with `include` to include all
      // rows which might otherwise have been filtered.
      val filters = optionalFilters.map(_.getOrElse(include))
      // AND together the filters.  If there are no filters at all return `include`.
      filters.reduceLeftOption(_ && _).getOrElse(include)
    }
  }

  def countWorkflowMetadataSummaryEntries(workflowStatuses: Set[String], workflowNames: Set[String],
                                          workflowExecutionUuids: Set[String], startTimestampOption: Option[Timestamp],
                                          endTimestampOption: Option[Timestamp]) = {
    val filter = filterWorkflowMetadataSummaryEntries(
      workflowStatuses, workflowNames, workflowExecutionUuids, startTimestampOption, endTimestampOption)
    workflowMetadataSummaryEntries.filter(filter).length
  }

  /**
    * Query workflow execution using the filter criteria encapsulated by the `WorkflowExecutionQueryParameters`.
    */
  def queryWorkflowMetadataSummaryEntries(workflowStatuses: Set[String], workflowNames: Set[String],
                                          workflowExecutionUuids: Set[String], startTimestampOption: Option[Timestamp],
                                          endTimestampOption: Option[Timestamp], page: Option[Int],
                                          pageSize: Option[Int]) = {
    val filter = filterWorkflowMetadataSummaryEntries(
      workflowStatuses, workflowNames, workflowExecutionUuids, startTimestampOption, endTimestampOption)
    val query = workflowMetadataSummaryEntries.filter(filter)
    (page, pageSize) match {
      case (Some(p), Some(ps)) => query.drop((p - 1) * ps).take(ps)
      case _ => query
    }
  }
}
