package cromwell.database.sql.tables

case class JobStoreResultSimpletonEntry
(
  simpletonKey: String,
  simpletonValue: String,
  wdlType: String,
  jobStoreId: Int,
  jobStoreSimpletonEntryId: Option[Int] = None
) extends DatabaseSimpleton

