package cromwell.database.slick

import java.sql.Timestamp

import cromwell.database.obj.FailureEvent

trait FailureEventComponent {
  this: DriverComponent with ExecutionComponent with WorkflowExecutionComponent =>

  import driver.api._

  class FailureEvents(tag: Tag) extends Table[FailureEvent](tag, "FAILURE_EVENT") {
    def failureId = column[Int]("FAILURE_EVENT_ID", O.PrimaryKey, O.AutoInc)
    def workflowExecutionId = column[Int]("WORKFLOW_EXECUTION_ID")
    def executionId = column[Option[Int]]("EXECUTION_ID")
    def message = column[String]("EVENT_MESSAGE")
    def timestamp = column[Timestamp]("EVENT_TIMESTAMP")

    override def * = (workflowExecutionId, executionId, message, timestamp, failureId.?) <>
      (FailureEvent.tupled, FailureEvent.unapply)

    def workflowExecution = foreignKey("FK_FAILURE_EVENT_WORKFLOW_EXECUTION_ID", workflowExecutionId, workflowExecutions)(_.workflowExecutionId)
    def execution = foreignKey("FK_FAILURE_EVENT_EXECUTION_ID", executionId, executions)(_.executionId.?)
  }

  protected val failureEvents = TableQuery[FailureEvents]

  val failureEventIdsAutoInc = failureEvents returning failureEvents.map(_.failureId)

  // Convenience function
  val failureEventsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
      (failureEvent, maybeExecution) <- failureEvents joinLeft executions on {
        case (failureEvent, execution) => failureEvent.executionId === execution.executionId
      }
      workflowExecution <- failureEvent.workflowExecution
      if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
    } yield (workflowExecution.workflowExecutionUuid, failureEvent.message, failureEvent.timestamp,
      maybeExecution.map(_.callFqn), maybeExecution.map(_.index), maybeExecution.map(_.attempt))
  )
}
