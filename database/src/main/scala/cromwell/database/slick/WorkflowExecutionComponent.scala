package cromwell.database.slick

import java.sql.Timestamp

import cromwell.database.obj.WorkflowExecution

trait WorkflowExecutionComponent {
  this: DriverComponent =>

  import driver.api._

  class WorkflowExecutions(tag: Tag) extends Table[WorkflowExecution](tag, "WORKFLOW_EXECUTION") {
    def workflowExecutionId = column[Int]("WORKFLOW_EXECUTION_ID", O.PrimaryKey, O.AutoInc)

    def workflowExecutionUuid = column[String]("WORKFLOW_EXECUTION_UUID")

    def name = column[String]("WORKFLOW_NAME")

    def status = column[String]("STATUS")

    def startDt = column[Timestamp]("START_DT")

    def endDt = column[Option[Timestamp]]("END_DT")

    override def * = (workflowExecutionUuid, name, status, startDt, endDt, workflowExecutionId.?) <>
      (WorkflowExecution.tupled, WorkflowExecution.unapply)

    def uniqueKey = index("UK_WE_WORKFLOW_EXECUTION_UUID",
      workflowExecutionUuid, unique = true)

    def statusIdx = index("STATUS_IDX", status, unique = false)
  }

  protected val workflowExecutions = TableQuery[WorkflowExecutions]

  val workflowExecutionIdsAutoInc = workflowExecutions returning workflowExecutions.map(_.workflowExecutionId)

  val workflowExecutionIdsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowExecution.workflowExecutionId)

  val workflowExecutionStatusByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowExecution.status)

  val workflowExecutionStatusEndDtByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (workflowExecution.status, workflowExecution.endDt))

  val workflowExecutionsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowExecution)

  /**
   * Query workflow execution using the filter criteria encapsulated by the `WorkflowExecutionQueryParameters`.
   */
  def queryWorkflowExecutions(statuses: Set[String], names: Set[String], startDate: Option[Timestamp],
                              endDate: Option[Timestamp]) = {
    val include: Rep[Boolean] = true
    val exclude: Rep[Boolean] = false
    workflowExecutions filter { workflow =>
      // All query parameters are either Options or Sets, so they might have no values specified at all.  The general
      // pattern for these criteria is to map Options and map/reduceLeftOption Sets, resulting in optional filters.

      // Start date is a non-null field and at most single-valued in the query parameters.
      val startDateTimeFilter = startDate.map(start => workflow.startDt >= start)
      // End date is nullable, necessitating the fold.  If the end date is null in the database we want to filter the
      // row if an end date filter has been specified.
      val endDateTimeFilter = endDate.map(end => workflow.endDt.fold(ifEmpty = exclude) { _ <= end })
      // Names and statuses are potentially multi-valued, the reduceLeftOption ORs together any name or status criteria
      // to include all matching names and statuses.
      val nameFilter = names.map(name => workflow.name === name).reduceLeftOption(_ || _)
      val statusFilter = statuses.map(status => workflow.status === status).reduceLeftOption(_ || _)

      // Put all the optional filters above together in one place.
      val optionalFilters: List[Option[Rep[Boolean]]] = List(nameFilter, statusFilter, startDateTimeFilter, endDateTimeFilter)
      // Unwrap the optional filters.  If any of these filters are not defined, replace with `include` to include all
      // rows which might otherwise have been filtered.
      val filters = optionalFilters.map(_.getOrElse(include))
      // AND together the filters.  If there are no filters at all return `include`.
      filters.reduceLeftOption(_ && _).getOrElse(include)
    }
  }
}
