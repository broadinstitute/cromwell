package cromwell.database.sql.tables

import java.sql.Clob

case class JobStoreEntry
(
  workflowExecutionUuid: String,
  callFullyQualifiedName: String,
  jobIndex: Int,
  jobAttempt: Int,
  jobSuccessful: Boolean,
  returnCode: Option[Int],
  exceptionMessage: Option[Clob],
  retryableFailure: Option[Boolean],
  jobStoreEntryId: Option[Int] = None
)
