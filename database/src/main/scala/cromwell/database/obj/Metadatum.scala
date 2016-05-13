package cromwell.database.obj

import java.sql.Timestamp

case class Metadatum
(
  workflowUuid: String,
  key: String,
  callFqn: Option[String],
  index: Option[Int],
  attempt: Option[Int],
  value: Option[String],
  timestamp: Timestamp,
  metadatumId: Option[Int] = None
)
