package cromwell.engine.db.slick

import java.sql.Timestamp

import cromwell.engine.ExecutionIndex._
import cromwell.engine.ExecutionStatus
import cromwell.engine.db.ExecutionDatabaseKey
import slick.profile.RelationalProfile.ColumnOption.Default

import scala.language.postfixOps

case class Execution(workflowExecutionId: Int,
                     callFqn: String,
                     index: Int,
                     status: String,
                     rc: Option[Int] = None,
                     startDt: Option[Timestamp] = None,
                     endDt: Option[Timestamp] = None,
                     backendType: String,
                     allowsResultReuse: Boolean = true,
                     dockerImageHash: Option[String] = None,
                     resultsClonedFrom: Option[Int] = None,
                     overallHash: Option[String] = None,
                     attempt: Int = 1,
                     executionId: Option[Int] = None)

trait ExecutionComponent {
  this: DriverComponent with WorkflowExecutionComponent =>

  import driver.api._

  class Executions(tag: Tag) extends Table[Execution](tag, "EXECUTION") {
    def executionId = column[Int]("EXECUTION_ID", O.PrimaryKey, O.AutoInc)
    def workflowExecutionId = column[Int]("WORKFLOW_EXECUTION_ID")
    def callFqn = column[String]("CALL_FQN")
    def index = column[Int]("IDX")
    def status = column[String]("STATUS")
    def rc = column[Option[Int]]("RC")
    def startDt = column[Option[Timestamp]]("START_DT")
    def endDt = column[Option[Timestamp]]("END_DT")
    def backendType = column[String]("BACKEND_TYPE")
    def allowsResultReuse = column[Boolean]("ALLOWS_RESULT_REUSE", Default(true))
    def dockerImageHash = column[Option[String]]("DOCKER_IMAGE_HASH")
    def resultsClonedFrom = column[Option[Int]]("RESULTS_CLONED_FROM")
    def executionHash = column[Option[String]]("EXECUTION_HASH")
    def attempt = column[Int]("ATTEMPT")

    override def * = (workflowExecutionId, callFqn, index, status, rc, startDt, endDt, backendType, allowsResultReuse, dockerImageHash, resultsClonedFrom, executionHash, attempt, executionId.?) <>
      (Execution.tupled, Execution.unapply)

    def workflowExecution = foreignKey(
      "FK_EXECUTION_WORKFLOW_EXECUTION_ID", workflowExecutionId, workflowExecutions)(_.workflowExecutionId)

    def resultsClonedFromExecution = foreignKey(
      "FK_RESULTS_CLONED_FROM", resultsClonedFrom, executions)(_.executionId.?, onDelete = ForeignKeyAction.SetNull)

    def hashIndex = index("HASH_INDEX", executionHash, unique = false)

    def uniqueKey = index("UK_WORKFLOW_CALL_INDEX_ATTEMPT", (workflowExecutionId, callFqn, index, attempt), unique = true)
  }

  protected val executions = TableQuery[Executions]

  val executionsAutoInc = executions returning executions.
    map(_.executionId) into ((a, id) => a.copy(executionId = Some(id)))

  val executionCallFqnsAndStatusesByWorkflowExecutionId = Compiled(
    (workflowExecutionId: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
    } yield execution)

  val executionStatusesAndReturnCodesByWorkflowExecutionIdAndCallKey = Compiled(
    (workflowExecutionId: Rep[Int], callFqn: Rep[String], index: Rep[Int], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
      if execution.callFqn === callFqn
      if execution.index === index
      if execution.attempt === attempt
    } yield (execution.status, execution.rc, execution.executionHash, execution.dockerImageHash))

  val executionStatusByWorkflowExecutionIdAndCallFqn = Compiled(
    (workflowExecutionId: Rep[Int], callFqn: Rep[String]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
      if execution.callFqn === callFqn
    } yield execution)

  val executionsByWorkflowExecutionUuidAndCallFqnAndShardIndexAndAttempt = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], index: Rep[Int], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.callFqn === callFqn
      if execution.index === index
      if execution.attempt === attempt
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield execution)

  val executionsByWorkflowExecutionIdAndCallFqnAndIndexAndAttempt = Compiled(
    (workflowExecutionId: Rep[Int], callFqn: Rep[String], index: Rep[Int], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
      if execution.callFqn === callFqn
      if execution.index === index
      if execution.attempt === attempt
    } yield execution)

  val executionsByWorkflowExecutionIdAndCallFqnAndAttempt = Compiled(
    (workflowExecutionId: Rep[Int], callFqn: Rep[String], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
      if execution.callFqn === callFqn
      if execution.attempt === attempt
    } yield execution)

  val executionsByWorkflowExecutionId = Compiled(
    (workflowExecutionId: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
    } yield execution)

  val executionsByWorkflowExecutionUuidAndCallFqnAndAttempt = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], attempt: Rep[Int]) => for {
      execution <- executions
      if execution.callFqn === callFqn
      if execution.attempt === attempt
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield execution)

  val executionsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      execution <- executions
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield execution)

  val executionsForRestartByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      execution <- executions
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
      if !(execution.status === ExecutionStatus.NotStarted.toString || execution.status === ExecutionStatus.Done.toString)
    } yield execution)

  val executionsWithReusableResultsByExecutionHash = Compiled(
    (executionHash: Rep[String]) => for {
      execution <- executions
      if execution.executionHash === executionHash && execution.allowsResultReuse
    } yield execution)

  val executionStatusesAndReturnCodesByExecutionId = Compiled(
    (executionId: Rep[Int]) => for {
      execution <- executions
      if execution.executionId === executionId
    } yield (execution.status, execution.rc))

  // see workflowExecutionsByStatuses
  def executionsByWorkflowExecutionIdAndScopeKeys(workflowExecutionId: Int, scopeKeys: Traversable[ExecutionDatabaseKey]) = {
    val falseRep: Rep[Boolean] = false
    val workflowFilteredQuery = for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
    } yield execution

    val scopeID = scopeKeys.map(k => (k.fqn, k.index.fromIndex, k.attempt)).toIterable

    /*
     * FIXME: This is bad, there is probably a better way
     * We want entries that have the right workflowID AND a (fqn, index) pair which is in scopeKeys
     * Use workaround because slick does not supporting tuples with inSet ATM: https://github.com/slick/slick/pull/995
     */
    workflowFilteredQuery filter { exec =>
      scopeID.map({
        case (name, index, attempt) => exec.callFqn === name && exec.index === index && exec.attempt === attempt
      }).fold(falseRep)(_ || _)
    }
  }
}
