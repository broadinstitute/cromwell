package cromwell.database.sql.tables

case class JobStoreEntry
(
  workflowExecutionUuid: String,
  callFullyQualifiedName: String,
  jobIndex: Int,
  jobAttempt: Int,
  jobSuccessful: Boolean,
  returnCode: Option[Int],
  exceptionMessage: Option[String],
  retryableFailure: Option[Boolean],
  jobStoreEntryId: Option[Int] = None
)
