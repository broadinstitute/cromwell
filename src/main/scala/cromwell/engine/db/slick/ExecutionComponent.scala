package cromwell.engine.db.slick

import java.sql.Timestamp

import cromwell.engine.ExecutionIndex._
import cromwell.engine.ExecutionStatus
import cromwell.engine.db.ExecutionDatabaseKey

import scala.language.postfixOps

case class Execution(workflowExecutionId: Int,
                     callFqn: String,
                     index: Int,
                     status: String,
                     rc: Option[Int] = None,
                     startDt: Option[Timestamp] = None,
                     endDt: Option[Timestamp] = None,
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

    override def * = (workflowExecutionId, callFqn, index, status, rc, startDt, endDt, executionId.?) <>
      (Execution.tupled, Execution.unapply)

    def workflowExecution = foreignKey(
      "FK_EXECUTION_WORKFLOW_EXECUTION_ID", workflowExecutionId, workflowExecutions)(_.workflowExecutionId)

    def uniqueKey = index("UK_WORKFLOW_CALL_INDEX", (workflowExecutionId, callFqn, index), unique = true)
  }

  protected val executions = TableQuery[Executions]

  val executionsAutoInc = executions returning executions.
    map(_.executionId) into ((a, id) => a.copy(executionId = Some(id)))

  val executionCallFqnsAndStatusesByWorkflowExecutionId = Compiled(
    (workflowExecutionId: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
    } yield (execution.callFqn, execution.index, execution.status, execution.rc))

  val executionStatusesAndReturnCodesByWorkflowExecutionIdAndCallKey = Compiled(
    (workflowExecutionId: Rep[Int], callFqn: Rep[String], index: Rep[Int]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
      if execution.callFqn === callFqn
      if execution.index === index
    } yield (execution.status, execution.rc))

  val executionStatusByWorkflowExecutionIdAndCallFqn = Compiled(
    (workflowExecutionId: Rep[Int], callFqn: Rep[String]) => for {
      execution <- executions
      if execution.workflowExecutionId === workflowExecutionId
      if execution.callFqn === callFqn
    } yield (execution.callFqn, execution.index, execution.status, execution.rc))

  val executionsByWorkflowExecutionUuidAndCallFqnAndShardIndex = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String], index: Rep[Int]) => for {
      execution <- executions
      if execution.callFqn === callFqn
      if execution.index === index
      workflowExecution <- execution.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield execution)

  val executionsByWorkflowExecutionUuidAndCallFqn = Compiled(
    (workflowExecutionUuid: Rep[String], callFqn: Rep[String]) => for {
      execution <- executions
      if execution.callFqn === callFqn
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

    val scopeID = scopeKeys.map(k => (k.fqn, k.index.fromIndex)).toIterable

    /*
     * FIXME: This is bad, there is probably a better way
     * We want entries that have the right workflowID AND a (fqn, index) pair which is in scopeKeys
     * Use workaround because slick does not supporting tuples with inSet ATM: https://github.com/slick/slick/pull/995 
     */
    workflowFilteredQuery filter { exec =>
      scopeID.map({
        case (name, index) => exec.callFqn === name && exec.index === index
      }).fold(falseRep)(_ || _)
    }
  }
}
