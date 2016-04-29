package cromwell.database.obj

import java.sql.Timestamp

case class ExecutionEvent
(
  executionId: Int,
  description: String,
  startTime: Timestamp,
  endTime: Timestamp,
  executionEventId: Option[Int] = None
)
