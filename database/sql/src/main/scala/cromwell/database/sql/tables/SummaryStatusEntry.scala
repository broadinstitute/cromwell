package cromwell.database.sql.tables

case class SummaryStatus
(
  summaryTableName: String,
  summarizedTableName: String,
  maximumId: Long,
  summaryStatusId: Option[Int] = None
)
