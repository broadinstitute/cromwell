package cromwell.database.sql.tables

case class CallCachingDetritusEntry
(
  detritusKey: String,
  detritusValue: String,
  callCachingEntryId: Option[Int] = None,
  callCachingDetritusId: Option[Int] = None
)
