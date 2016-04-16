package cromwell.database.obj

import java.sql.Clob

case class Symbol
(
  workflowExecutionId: Int,
  scope: String,
  name: String,
  index: Int, // https://bugs.mysql.com/bug.php?id=8173
  io: String,
  reportableResult: Boolean,
  wdlType: String,
  wdlValue: Option[Clob],
  symbolHash: Option[String],
  symbolId: Option[Int] = None
)
