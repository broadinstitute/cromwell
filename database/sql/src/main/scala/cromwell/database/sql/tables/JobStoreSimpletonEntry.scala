package cromwell.database.sql.tables

import javax.sql.rowset.serial.SerialClob

case class JobStoreSimpletonEntry
(
  simpletonKey: String,
  simpletonValue: Option[SerialClob],
  wdlType: String,
  jobStoreEntryId: Option[Long] = None,
  jobStoreSimpletonEntryId: Option[Long] = None
)
