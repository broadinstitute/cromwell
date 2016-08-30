package cromwell.database.sql.tables

case class CallCachingHashEntry
(
  hashKey: String,
  hashValue: String,
  resultMetaInfoId: Int,
  callCachingHashId: Option[Int] = None
)
