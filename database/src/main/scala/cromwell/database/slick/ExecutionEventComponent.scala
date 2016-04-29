package cromwell.database.slick

import java.sql.Timestamp

import cromwell.database.obj.ExecutionEvent

trait ExecutionEventComponent {
  this: DriverComponent with ExecutionComponent =>

  import driver.api._

  class ExecutionEvents(tag: Tag) extends Table[ExecutionEvent](tag, "EXECUTION_EVENT") {
    def executionEventId = column[Int]("EVENT_ID", O.PrimaryKey, O.AutoInc)
    def executionId = column[Int]("EXECUTION_ID")
    def description = column[String]("DESCRIPTION")
    def startTime = column[Timestamp]("START_DT")
    def endTime = column[Timestamp]("END_DT")

    override def * = (executionId, description, startTime, endTime, executionEventId.?) <>
      (ExecutionEvent.tupled, ExecutionEvent.unapply)

    def execution = foreignKey("FK_EXECUTION_EVENT_EXECUTION_ID", executionId, executions)(_.executionId)
    def uniqueKey = index("UK_EXECUTION_EVENT_EXECUTION_ID_DESCRIPTION", (executionId, description), unique = true)
  }

  protected val executionEvents = TableQuery[ExecutionEvents]

  val executionEventIdsAutoInc = executionEvents returning executionEvents.map(_.executionEventId)

  // Convenience function
  val executionEventsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
        executionEvent <- executionEvents
        execution <- executionEvent.execution
        workflowExecution <- execution.workflowExecution
        if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
      } yield (execution.callFqn, execution.index, execution.attempt,
      executionEvent.description, executionEvent.startTime, executionEvent.endTime)
  )
}
