package cromwell.database.sql.tables

case class SummaryStatusEntry
(
  summaryName: String,
  summaryPosition: Long,
  summaryStatusEntryId: Option[Int] = None
)
