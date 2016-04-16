package cromwell.database.obj

import java.sql.Timestamp

case class FailureEvent
(
  workflowExecutionId: Int,
  executionId: Option[Int],
  failure: String,
  timestamp: Timestamp,
  failureId: Option[Int] = None
)
