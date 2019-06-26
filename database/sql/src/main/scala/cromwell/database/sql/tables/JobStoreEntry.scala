package cromwell.database.sql.tables

import javax.sql.rowset.serial.SerialClob

case class JobStoreEntry
(
  workflowExecutionUuid: String,
  callFullyQualifiedName: String,
  jobIndex: Int,
  jobAttempt: Int,
  jobSuccessful: Boolean,
  returnCode: Option[Int],
  exceptionMessage: Option[SerialClob],
  retryableFailure: Option[Boolean],
  jobStoreEntryId: Option[Int] = None
)
