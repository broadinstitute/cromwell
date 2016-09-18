package cromwell.database.slick.tables

import java.sql.Timestamp

import cromwell.database.sql.tables.WorkflowMetadataSummary

trait WorkflowMetadataSummaryComponent {

  this: DriverComponent with MetadataComponent =>

  import driver.api._

  class WorkflowMetadataSummaries(tag: Tag) extends Table[WorkflowMetadataSummary](tag, "WORKFLOW_METADATA_SUMMARY") {
    def workflowMetadataSummaryId = column[Long]("WORKFLOW_METADATA_SUMMARY_ID", O.PrimaryKey, O.AutoInc)
    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")
    def name = column[Option[String]]("WORKFLOW_NAME")
    def status = column[Option[String]]("WORKFLOW_STATUS")
    def startDate = column[Option[Timestamp]]("START_DT")
    def endDate = column[Option[Timestamp]]("END_DT")

    // WorkflowMetadataSummary has a companion object which apparently necessitates the weird apply / tupled syntax.
    override def * = (workflowExecutionUuid, name, status, startDate, endDate, workflowMetadataSummaryId.?) <>
      ((WorkflowMetadataSummary.apply _ ).tupled, WorkflowMetadataSummary.unapply)

    def uuidIndex = index("WORKFLOW_METADATA_UUID_IDX", workflowExecutionUuid, unique = true)

    def nameIndex = index("WORKFLOW_METADATA_NAME_IDX", name, unique = false)

    def statusIndex = index("WORKFLOW_METADATA_STATUS_IDX", status, unique = false)
  }

  val workflowMetadataSummaries = TableQuery[WorkflowMetadataSummaries]

  val workflowMetadataSummaryAutoInc = workflowMetadataSummaries returning workflowMetadataSummaries.map(_.workflowMetadataSummaryId)

  val workflowMetadataSummariesByUuid = Compiled(
    (workflowUuid: Rep[String]) => for {
      workflowMetadataSummary <- workflowMetadataSummaries
      if workflowMetadataSummary.workflowExecutionUuid === workflowUuid
    } yield workflowMetadataSummary)

  val workflowStatusByUuid = Compiled(
    (workflowUuid: Rep[String]) => for {
      workflowMetadataSummary <- workflowMetadataSummaries
      if workflowMetadataSummary.workflowExecutionUuid === workflowUuid
    } yield workflowMetadataSummary.status
  )

  def filterWorkflowSummaries(statuses: Set[String], names: Set[String], uuids: Set[String],
                              startDate: Option[Timestamp], endDate: Option[Timestamp]): (WorkflowMetadataSummaries) => Rep[Boolean] = {
    val include: Rep[Boolean] = true
    val exclude: Rep[Boolean] = false

    { workflowSummary =>
      // All query parameters are either Options or Sets, so they might have no values specified at all.  The general
      // pattern for these criteria is to map Options and map/reduceLeftOption Sets, resulting in optional filters.

      // All fields but UUID are nullable, necessitating the folds.  If these fields are null in the database we want to
      // filter the row if the relevant filter has been specified.
      val startDateTimeFilter = startDate.map(start => workflowSummary.startDate.fold(ifEmpty = exclude) { _ >= start })
      val endDateTimeFilter = endDate.map(end => workflowSummary.endDate.fold(ifEmpty = exclude) { _ <= end })
      // Names, UUIDs, and statuses are potentially multi-valued, the reduceLeftOption ORs together any name, UUID, or
      // status criteria to include all matching names, UUIDs, and statuses.
      val nameFilter = names.map(name => workflowSummary.name.fold(ifEmpty = exclude) { _ === name }).reduceLeftOption(_ || _)
      val uuidFilter = uuids.map(uuid => workflowSummary.workflowExecutionUuid === uuid).reduceLeftOption(_ || _)
      val statusFilter = statuses.map(status => workflowSummary.status.fold(ifEmpty = exclude) { _ === status }).reduceLeftOption(_ || _)

      // Put all the optional filters above together in one place.
      val optionalFilters: List[Option[Rep[Boolean]]] =
        List(nameFilter, uuidFilter, statusFilter, startDateTimeFilter, endDateTimeFilter)
      // Unwrap the optional filters.  If any of these filters are not defined, replace with `include` to include all
      // rows which might otherwise have been filtered.
      val filters = optionalFilters.map(_.getOrElse(include))
      // AND together the filters.  If there are no filters at all return `include`.
      filters.reduceLeftOption(_ && _).getOrElse(include)
    }
  }

  def countWorkflowSummaries(statuses: Set[String], names: Set[String], uuids: Set[String],
                             startDate: Option[Timestamp], endDate: Option[Timestamp]) = {
    workflowMetadataSummaries.filter(filterWorkflowSummaries(statuses, names, uuids, startDate, endDate)).length
  }

  /**
    * Query workflow execution using the filter criteria encapsulated by the `WorkflowExecutionQueryParameters`.
    */
  def queryWorkflowSummaries(statuses: Set[String], names: Set[String], uuids: Set[String],
                             startDate: Option[Timestamp], endDate: Option[Timestamp],
                             page: Option[Int], pageSize: Option[Int]) = {
    val query = workflowMetadataSummaries.filter(filterWorkflowSummaries(statuses, names, uuids, startDate, endDate))
    (page, pageSize) match {
      case (Some(p), Some(ps)) => query.drop((p - 1) * ps).take(ps)
      case _ => query
    }
  }
}
