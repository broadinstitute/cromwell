package cromwell.database.sql.tables

import java.sql.Timestamp

import javax.sql.rowset.serial.SerialClob

case class MetadataJournalEntry
(
  workflowExecutionUuid: String,
  callFullyQualifiedName: Option[String],
  jobIndex: Option[Int],
  jobAttempt: Option[Int],
  metadataKey: String,
  metadataValue: Option[SerialClob],
  metadataValueType: Option[String],
  metadataGenerationTimestamp: Timestamp,
  metadataWriteTimestamp: Timestamp,
  summarizationTimestamp: Option[Timestamp],
  metadataEntryId: Option[Long] = None
)
