package cromwell.database.sql.tables

case class CallCachingSimpletonEntry
(
  simpletonKey: String,
  simpletonValue: String,
  wdlType: String,
  callCachingEntryId: Option[Int] = None,
  callCachingSimpletonEntryId: Option[Int] = None
)
