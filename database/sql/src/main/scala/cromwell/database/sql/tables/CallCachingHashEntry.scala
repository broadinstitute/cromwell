package cromwell.database.sql.tables

case class CallCachingHashEntry
(
  hashKey: String,
  hashValue: String,
  callCachingEntryId: Option[Int] = None,
  callCachingHashEntryId: Option[Int] = None
)
