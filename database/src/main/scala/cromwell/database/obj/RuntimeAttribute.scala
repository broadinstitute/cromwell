package cromwell.database.obj

case class RuntimeAttribute
(
  executionId: Int,
  name: String,
  value: String,
  runtimeAttributeId: Option[Int] = None
)
