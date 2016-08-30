package cromwell.database.sql.tables

import java.sql.Timestamp

@deprecated("Olde Worlde Databasee Tablee", "0.21")
case class WorkflowExecution
(
  workflowExecutionUuid: String,
  name: String,
  status: String,
  startDt: Timestamp,
  endDt: Option[Timestamp],
  workflowExecutionId: Option[Int] = None
)
