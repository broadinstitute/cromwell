package cromwell.database.sql.tables

import java.sql.Timestamp

import javax.sql.rowset.serial.SerialClob

case class MetadataJournalEntry
(
  workflowExecutionUuid: String,
  callFqn: Option[String],
  jobScatterIndex: Option[Int],
  jobRetryAttempt: Option[Int],
  metadataKey: String,
  metadataValue: Option[SerialClob],
  metadataValueType: Option[String],
  metadataTimestamp: Timestamp,
  needsSummarization: Boolean,
  metadataJournalId: Option[Long] = None
) {
  def asMetadataEntry = MetadataEntry(
    workflowExecutionUuid = workflowExecutionUuid,
    callFullyQualifiedName = callFqn,
    jobIndex = jobScatterIndex,
    jobAttempt = jobRetryAttempt,
    metadataKey = metadataKey,
    metadataValue = metadataValue,
    metadataValueType = metadataValueType,
    metadataTimestamp = metadataTimestamp
  )
}
