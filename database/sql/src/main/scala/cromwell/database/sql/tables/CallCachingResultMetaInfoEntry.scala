package cromwell.database.sql.tables

case class CallCachingResultMetaInfoEntry
(
  workflowUuid: String,
  callFqn: String,
  scatterIndex: Int,
  returnCode: Option[Int],
  allowResultReuse: Boolean,
  callCachingResultMetaInfoEntryId: Option[Int] = None
)
