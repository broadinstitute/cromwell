package cromwell.database.sql.tables

import java.sql.{Clob, Timestamp}

case class MetadataEntry
(
  workflowExecutionUuid: String,
  callFullyQualifiedName: Option[String],
  jobIndex: Option[Int],
  jobAttempt: Option[Int],
  metadataKey: String,
  metadataValue: Option[Clob],
  metadataValueType: Option[String],
  metadataTimestamp: Timestamp,
  metadataEntryId: Option[Long] = None
)
