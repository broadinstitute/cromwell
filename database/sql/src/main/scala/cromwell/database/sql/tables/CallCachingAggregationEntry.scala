package cromwell.database.sql.tables

case class CallCachingAggregationEntry
(
  baseAggregation: String,
  inputFilesAggregation: Option[String],
  callCachingEntryId: Option[Long] = None,
  callCachingAggregationEntryId: Option[Int] = None
)
