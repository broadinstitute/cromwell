package cromwell.database.slick.tables

import java.sql.Timestamp

import cromwell.database.sql.tables.Execution
import slick.jdbc.GetResult
import slick.profile.RelationalProfile.ColumnOption.Default

@deprecated("Olde Worlde Databasee Tablee", "0.21")
trait ExecutionComponent {
  this: DriverComponent with WorkflowExecutionComponent =>

  import driver.api._

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  class Executions(tag: Tag) extends Table[Execution](tag, "EXECUTION") {
    def executionId = column[Int]("EXECUTION_ID", O.PrimaryKey, O.AutoInc)
    def workflowExecutionId = column[Int]("WORKFLOW_EXECUTION_ID")
    def callFqn = column[String]("CALL_FQN")
    def index = column[Int]("IDX")
    def attempt = column[Int]("ATTEMPT")
    def status = column[String]("STATUS")
    def rc = column[Option[Int]]("RC")
    def startDt = column[Option[Timestamp]]("START_DT")
    def endDt = column[Option[Timestamp]]("END_DT")
    def backendType = column[String]("BACKEND_TYPE")
    def allowsResultReuse = column[Boolean]("ALLOWS_RESULT_REUSE", Default(true))
    def dockerImageHash = column[Option[String]]("DOCKER_IMAGE_HASH")
    def resultsClonedFrom = column[Option[Int]]("RESULTS_CLONED_FROM")
    def executionHash = column[Option[String]]("EXECUTION_HASH")

    override def * = (workflowExecutionId, callFqn, index, attempt, status, rc, startDt, endDt, backendType,
      allowsResultReuse, dockerImageHash, resultsClonedFrom, executionHash, executionId.?) <>
      (Execution.tupled, Execution.unapply)

    def workflowExecution = foreignKey(
      "FK_EXECUTION_WORKFLOW_EXECUTION_ID", workflowExecutionId, workflowExecutions)(_.workflowExecutionId)

    def resultsClonedFromExecution = foreignKey(
      "FK_RESULTS_CLONED_FROM", resultsClonedFrom, executions)(_.executionId.?, onDelete = ForeignKeyAction.SetNull)

    def hashIndex = index("HASH_INDEX", executionHash, unique = false)

    def uniqueKey = index("UK_WORKFLOW_CALL_INDEX_ATTEMPT", (workflowExecutionId, callFqn, index, attempt), unique = true)
  }

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  protected val executions = TableQuery[Executions]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionIdsAutoInc = executions returning executions.map(_.executionId)

  private[this] implicit val executionResult = GetResult { r =>
    Execution(r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<, r.<<)
  }

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionCallStatusesByWorkflowExecutionUuidAndCallKey = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], index: Rep[Int], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.callFqn === callFqn
      if execution.index === index
      if execution.attempt === attempt
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (execution.status, execution.rc, execution.executionHash, execution.dockerImageHash))

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionCallStatusesByWorkflowExecutionUuidAndCallFqn = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String]) => for {
      execution <- executions
      if execution.callFqn === callFqn
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (execution.callFqn, execution.index, execution.attempt,
      execution.status, execution.rc, execution.executionHash, execution.dockerImageHash))

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionCallStatusesByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      execution <- executions
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (execution.callFqn, execution.index, execution.attempt,
      execution.status, execution.rc, execution.executionHash, execution.dockerImageHash))

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionIdsByWorkflowExecutionUuidAndCallKey = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], index: Rep[Int], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.callFqn === callFqn
      if execution.index === index
      if execution.attempt === attempt
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield execution.executionId)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionIdAndWorkflowExecutionIdsByWorkflowExecutionUuidAndCallKey = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], index: Rep[Int], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.callFqn === callFqn
      if execution.index === index
      if execution.attempt === attempt
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (execution.executionId, workflowExecution.workflowExecutionId))

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionAllowsResultReusesByWorkflowExecutionIdAndCallKey = Compiled(
    (workflowExecutionId: Rep[Int], callFqn: Rep[String], index: Rep[Int], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
      if execution.callFqn === callFqn
      if execution.index === index
      if execution.attempt === attempt
    } yield execution.allowsResultReuse)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionAllowsResultReusesByWorkflowExecutionIdAndCallFqnAndAttempt = Compiled(
    (workflowExecutionId: Rep[Int], callFqn: Rep[String], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
      if execution.callFqn === callFqn
      if execution.attempt === attempt
    } yield execution.allowsResultReuse)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionAllowsResultReusesByWorkflowExecutionId = Compiled(
    (workflowExecutionId: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
    } yield execution.allowsResultReuse)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionIdsByWorkflowExecutionUuidAndCallFqnAndAttempt = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.callFqn === callFqn
      if execution.attempt === attempt
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield execution.executionId)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionIdAndWorkflowExecutionIdsByWorkflowExecutionUuidAndCallFqnAndAttempt = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.callFqn === callFqn
      if execution.attempt === attempt
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (execution.executionId, workflowExecution.workflowExecutionId))

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      execution <- executions
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield execution)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionsWithReusableResultsByExecutionHash = Compiled(
    (executionHash: Rep[String]) => for {
      execution <- executions
      if execution.executionHash === executionHash
      if execution.allowsResultReuse
    } yield execution)

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def selectExecutionsWithCacheHitWorkflowAndCall(workflowExecutionUuid: String):
  DBIO[Vector[(Execution, Option[String], Option[String])]] = {
    /** NOTE: the columns must be in the same order as the attributes in the Execution class for this to work */
    sql"""
       SELECT e.WORKFLOW_EXECUTION_ID, e.CALL_FQN, e.IDX, e.ATTEMPT, e.STATUS, e.RC, e.START_DT, e.END_DT,
              e.BACKEND_TYPE, e.ALLOWS_RESULT_REUSE, e.DOCKER_IMAGE_HASH,
              e.RESULTS_CLONED_FROM, e.EXECUTION_HASH, e.EXECUTION_ID,
              we1.WORKFLOW_EXECUTION_UUID, e1.CALL_FQN
       FROM EXECUTION e
       LEFT JOIN EXECUTION e1 ON e.RESULTS_CLONED_FROM=e1.EXECUTION_ID
       LEFT JOIN WORKFLOW_EXECUTION we1 ON e1.WORKFLOW_EXECUTION_ID=we1.WORKFLOW_EXECUTION_ID
       LEFT JOIN WORKFLOW_EXECUTION we ON e.WORKFLOW_EXECUTION_ID=we.WORKFLOW_EXECUTION_ID
       WHERE we.WORKFLOW_EXECUTION_UUID=$workflowExecutionUuid
       """.as[(Execution, Option[String], Option[String])]
  }

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionsByWorkflowExecutionIdAndCallKey = Compiled(
    (workflowExecutionId: Rep[Int], callFqn: Rep[String], index: Rep[Int], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
      if execution.callFqn === callFqn
      if execution.index === index
      if execution.attempt === attempt
    } yield (execution.status, execution.endDt, execution.rc,
      execution.executionHash, execution.dockerImageHash, execution.resultsClonedFrom)
  )

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionStatusesByWorkflowExecutionIdAndCallKey = Compiled(
    (workflowExecutionId: Rep[Int], callFqn: Rep[String], index: Rep[Int], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
      if execution.callFqn === callFqn
      if execution.index === index
      if execution.attempt === attempt
    } yield execution.status
  )

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  val executionStatusesAndStartDtByWorkflowExecutionIdAndCallKey = Compiled(
    (workflowExecutionId: Rep[Int], callFqn: Rep[String], index: Rep[Int], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
      if execution.callFqn === callFqn
      if execution.index === index
      if execution.attempt === attempt
    } yield (execution.status, execution.startDt)
  )
}
