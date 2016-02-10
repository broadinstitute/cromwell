package cromwell.engine.db.slick

import java.sql.{Timestamp, Clob}

import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.{WorkflowId, ExecutionIndex}
import cromwell.engine.ExecutionIndex._

case class ExecutionEvent (
  executionId: Int,
  description: String,
  startTime: Timestamp,
  endTime: Timestamp,
  executionEventId: Option[Int] = None
)

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

  val executionEventsAutoInc = executionEvents returning executionEvents.
    map(_.executionEventId) into ((a, id) => a.copy(executionEventId = Some(id)))

  // Convenience function
  def executionEventsByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
        executionEvent <- executionEvents
        execution <- executionEvent.execution
        workflowExecution <- execution.workflowExecution
        if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
      } yield ((execution.callFqn, execution.index, execution.attempt), executionEvent)
  )
}
