package cromwell.database.sql.tables

case class CallCachingHashEntry
(
  hashKey: String,
  hashValue: String,
  callCachingEntryId: Option[Long] = None,
  callCachingHashEntryId: Option[Long] = None
)
