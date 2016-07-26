package cromwell.database.sql.tables

case class JobStoreEntry
(
  workflowUuid: String,
  callFqn: String,
  index: Int,
  attempt: Int,
  jobSuccessful: Boolean,
  returnCode: Option[Int],
  jobResult: Option[String],
  exceptionMessage: Option[String],
  jobStoreId: Option[Int] = None
)
