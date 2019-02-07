package cromwell.database.sql.tables

case class DockerHashStoreEntry
(
  workflowExecutionUuid: String,
  dockerTag: String,
  dockerHash: String,
  dockerSize: Option[Long],
  dockerHashStoreEntryId: Option[Int] = None
)
