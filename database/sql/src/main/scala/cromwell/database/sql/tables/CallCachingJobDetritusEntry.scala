package cromwell.database.sql.tables

final case class CallCachingJobDetritusEntry
(
  jobDetritusKey: String,
  jobDetritusValue: String,
  resultMetaInfoId: Int,
  callCachingJobDetritusId: Option[Int] = None
)
