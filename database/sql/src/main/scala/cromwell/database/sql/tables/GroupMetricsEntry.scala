package cromwell.database.sql.tables

import java.sql.Timestamp

case class GroupMetricsEntry(
  groupId: String,
  quotaExhaustionDetected: Timestamp,
  groupMetricsEntryId: Option[Long] = None
)
