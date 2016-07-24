package cromwell.database.sql.tables

@deprecated("Olde Worlde Databasee Tablee", "0.21")
case class RuntimeAttribute
(
  executionId: Int,
  name: String,
  value: String,
  runtimeAttributeId: Option[Int] = None
)
