package cromwell.database.sql.tables

case class JobStoreSimpletonEntry
(
  simpletonKey: String,
  simpletonValue: String,
  wdlType: String,
  jobStoreEntryId: Option[Int] = None,
  jobStoreSimpletonEntryId: Option[Int] = None
)
