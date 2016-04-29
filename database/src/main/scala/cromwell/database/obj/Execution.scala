package cromwell.database.obj

import java.sql.Timestamp

case class Execution
(
  workflowExecutionId: Int,
  callFqn: String,
  index: Int,
  attempt: Int,
  status: String,
  rc: Option[Int],
  startDt: Option[Timestamp],
  endDt: Option[Timestamp],
  backendType: String,
  allowsResultReuse: Boolean,
  dockerImageHash: Option[String],
  resultsClonedFrom: Option[Int],
  overallHash: Option[String],
  executionId: Option[Int] = None
)
