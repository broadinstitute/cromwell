package cromwell.database.obj

case class ExecutionInfo
(
  executionId: Int,
  key: String,
  value: Option[String],
  executionInfoId: Option[Int] = None
)
