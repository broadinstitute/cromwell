package cromwell.database.sql.tables

case class JobStoreEntry
(
  workflowUuid: String,
  callFqn: String,
  index: Int,
  attempt: Int,
  jobSuccessful: Boolean,
  returnCode: Option[Int],
  exceptionMessage: Option[String],
  retryableFailure: Option[Boolean],
  jobStoreId: Option[Int] = None
)
