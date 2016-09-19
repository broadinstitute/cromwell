package cromwell.database.sql.tables

import java.sql.Timestamp

case class Metadatum
(
  workflowUuid: String,
  key: String,
  callFqn: Option[String],
  index: Option[Int],
  attempt: Option[Int],
  value: Option[String],
  valueType: Option[String],
  timestamp: Timestamp,
  metadatumId: Option[Long] = None
)
