package cromwell.database.sql.tables

case class DockerHashStoreEntry
(
  workflowExecutionUuid: String,
  dockerTag: String,
  dockerHash: String,
  dockerHashStoreEntryId: Option[Int] = None
)
