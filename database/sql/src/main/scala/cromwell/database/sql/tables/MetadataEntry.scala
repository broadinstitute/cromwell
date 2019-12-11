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
) {
  def asMetadataJournalEntry = MetadataJournalEntry(
    workflowExecutionUuid = workflowExecutionUuid,
    callFqn = callFullyQualifiedName,
    jobScatterIndex = jobIndex,
    jobRetryAttempt = jobAttempt,
    metadataKey = metadataKey,
    metadataValue = metadataValue,
    metadataValueType = metadataValueType,
    metadataTimestamp = metadataTimestamp,
    needsSummarization = true // TODO: hmm, could we be smarter about which keys start off needing summarization?
  )
}
