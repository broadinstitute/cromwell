package cromwell.database.sql.tables

import java.sql.Timestamp

case class SummaryQueueEntry
(
  metadataEntryId: Long,
  metadateEntryWriteTimestamp: Timestamp,
  summarizationTimestamp: Option[Timestamp],
  summaryQueueId: Option[Long] = None
)
