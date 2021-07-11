package cromwell.database.sql.tables

final case class InformationSchemaEntry
(
  dataLength: Long,
  indexLength: Long,
  dataFree: Long
)
