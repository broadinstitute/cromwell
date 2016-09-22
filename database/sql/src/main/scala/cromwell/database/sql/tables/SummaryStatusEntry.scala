package cromwell.database.sql.tables

case class SummaryStatusEntry
(
  summaryTableName: String,
  summarizedTableName: String,
  maximumId: Long,
  summaryStatusEntryId: Option[Int] = None
)
