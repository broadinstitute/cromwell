package cromwell.database.sql.tables

import java.sql.Timestamp
import java.time.Instant

case class CallCachingEntry(
  workflowExecutionUuid: String,
  callFullyQualifiedName: String,
  jobIndex: Int,
  jobAttempt: Option[Int],
  returnCode: Option[Int],
  allowResultReuse: Boolean,
  callCachingEntryId: Option[Long] = None,
  // timestamp is set to time of CallCachingEntry creation by default
  createdAt: Timestamp = Timestamp.from(Instant.now())
)
