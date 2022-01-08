package cromwell.database.sql.tables

import javax.sql.rowset.serial.SerialClob

case class CallCachingSimpletonEntry
(
  simpletonKey: String,
  simpletonValue: Option[SerialClob],
  wdlType: String,
  callCachingEntryId: Option[Int] = None,
  callCachingSimpletonEntryId: Option[Int] = None
)
