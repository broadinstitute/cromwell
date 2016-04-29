package cromwell.database.obj

import java.sql.Timestamp

case class WorkflowExecution
(
  workflowExecutionUuid: String,
  name: String,
  status: String,
  startDt: Timestamp,
  endDt: Option[Timestamp],
  workflowExecutionId: Option[Int] = None
)
