package cromwell.database.sql.tables

@deprecated("Olde Worlde Databasee Tablee", "0.21")
case class ExecutionInfo
(
  executionId: Int,
  key: String,
  value: Option[String],
  executionInfoId: Option[Int] = None
)
