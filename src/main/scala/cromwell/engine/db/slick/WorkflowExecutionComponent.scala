package cromwell.engine.db.slick

import java.sql.Timestamp

import cromwell.webservice.WorkflowQueryParameters

case class WorkflowExecution
(
  workflowExecutionUuid: String,
  name: String,
  status: String,
  startDt: Timestamp,
  workflowExecutionId: Option[Int] = None,
  endDt: Option[Timestamp] = None
  )

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

    override def * = (workflowExecutionUuid, name, status, startDt, workflowExecutionId.?, endDt) <>
      (WorkflowExecution.tupled, WorkflowExecution.unapply)

    def uniqueKey = index("UK_WE_WORKFLOW_EXECUTION_UUID",
      workflowExecutionUuid, unique = true)
  }

  protected val workflowExecutions = TableQuery[WorkflowExecutions]

  val workflowExecutionsAutoInc = workflowExecutions returning workflowExecutions.
    map(_.workflowExecutionId) into ((a, id) => a.copy(workflowExecutionId = Some(id)))

  val workflowExecutionsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      workflowExecution <- workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield workflowExecution)

  // NOTE: No precompile for you!
  // [Compile] works for all functions ... consisting only of individual columns
  // - http://slick.typesafe.com/doc/3.0.0/queries.html
  // - https://groups.google.com/forum/#!topic/scalaquery/2d_r4DEthfY
  //
  // On one hand:
  // a) This should be compilable if we narrowed the API to only get Submitted/Running workflows... BUT
  // b) It turns out this is currently executed only once at server restart
  def workflowExecutionsByStatuses(statuses: Traversable[String]) = {
    for {
      workflowExecution <- workflowExecutions
      if workflowExecution.status inSet statuses
    } yield workflowExecution
  }

  /**
   * Query workflow execution using the filter criteria encapsulated by the `WorkflowQueryParameters`.
   */
  def queryWorkflowExecutions(queryParameters: WorkflowQueryParameters) = {
    val include: Rep[Boolean] = true
    val exclude: Rep[Boolean] = false
    workflowExecutions.filter { workflow =>
      // All query parameters are either Options or Sets, so they might have no values specified at all.  The general
      // pattern for these criteria is to map Options and map/reduceLeftOption Sets, resulting in optional filters.

      // Start date is a non-null field and at most single-valued in the query parameters.
      val startDateTimeFilter = queryParameters.startDate.map(start => workflow.startDt >= new Timestamp(start.getMillis))
      // End date is nullable, necessitating the fold.  If the end date is null in the database we want to filter the
      // row if an end date filter has been specified.
      val endDateTimeFilter = queryParameters.endDate.map(end => workflow.endDt.fold(ifEmpty = exclude) { _ <= new Timestamp(end.getMillis) })
      // Names and statuses are potentially multi-valued, the reduceLeftOption ORs together any name or status criteria
      // to include all matching names and statuses.
      val nameFilter = queryParameters.names.map(name => workflow.name === name).reduceLeftOption(_ || _)
      val statusFilter = queryParameters.statuses.map(status => workflow.status === status).reduceLeftOption(_ || _)

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
