package cromwell.database.sql.tables

import javax.sql.rowset.serial.SerialClob

case class JobStoreSimpletonEntry
(
  simpletonKey: String,
  simpletonValue: Option[SerialClob],
  wdlType: String,
  jobStoreEntryId: Option[Int] = None,
  jobStoreSimpletonEntryId: Option[Int] = None
)
