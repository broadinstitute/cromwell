package cromwell.database.sql.tables

import java.sql.Timestamp

@deprecated("Olde Worlde Databasee Tablee", "0.21")
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
