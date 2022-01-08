package cromwell.database.sql.tables

import java.sql.Timestamp

import javax.sql.rowset.serial.SerialClob

case class MetadataEntry
(
  workflowExecutionUuid: String,
  callFullyQualifiedName: Option[String],
  jobIndex: Option[Int],
  jobAttempt: Option[Int],
  metadataKey: String,
  metadataValue: Option[SerialClob],
  metadataValueType: Option[String],
  metadataTimestamp: Timestamp,
  metadataEntryId: Option[Long] = None
)
