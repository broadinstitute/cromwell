package cromwell.database.sql.tables

case class CallCachingResultSimpletonEntry
(
  simpletonKey: String,
  simpletonValue: String,
  wdlType: String,
  resultMetaInfoId: Int,
  callCachingResultSimpletonId: Option[Int]
) extends DatabaseSimpleton
