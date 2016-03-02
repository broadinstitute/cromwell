package cromwell.engine.db.slick

import java.sql.Timestamp

case class FailureEvent (
  workflowExecutionId: Int,
  executionId: Option[Int],
  failure: String,
  timestamp: Timestamp,
  failureId: Option[Int] = None
)

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

  val failureEventsAutoInc = failureEvents returning failureEvents.
    map(_.failureId) into ((a, id) => a.copy(failureId = Option(id)))

  // Convenience function
  def failuresByWorkflowExecutionUuid = Compiled(
    (workflowExecutionUuid: Rep[String]) => for {
        failureEvent <- failureEvents
        workflowExecution <- failureEvent.workflowExecution
        if workflowExecution.workflowExecutionUuid === workflowExecutionUuid
      } yield (workflowExecution.workflowExecutionUuid, failureEvent)
  )
}