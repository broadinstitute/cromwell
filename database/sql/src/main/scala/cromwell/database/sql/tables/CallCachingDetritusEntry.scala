package cromwell.database.sql.tables

import javax.sql.rowset.serial.SerialClob

case class CallCachingDetritusEntry
(
  detritusKey: String,
  detritusValue: Option[SerialClob],
  callCachingEntryId: Option[Int] = None,
  callCachingDetritusEntryId: Option[Int] = None
)
