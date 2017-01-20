package cromwell.database.sql.tables

import java.sql.Clob

case class CallCachingDetritusEntry
(
  detritusKey: String,
  detritusValue: Clob,
  callCachingEntryId: Option[Int] = None,
  callCachingDetritusId: Option[Int] = None
)
