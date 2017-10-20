package cromwell.database.sql.tables

import java.sql.Clob

case class CallCachingSimpletonEntry
(
  simpletonKey: String,
  simpletonValue: Option[Clob],
  wdlType: String,
  callCachingEntryId: Option[Int] = None,
  callCachingSimpletonEntryId: Option[Int] = None
)
