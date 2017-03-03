package cromwell.database.sql.tables

import java.sql.Clob

case class JobStoreSimpletonEntry
(
  simpletonKey: String,
  simpletonValue: Option[Clob],
  wdlType: String,
  jobStoreEntryId: Option[Int] = None,
  jobStoreSimpletonEntryId: Option[Int] = None
)
