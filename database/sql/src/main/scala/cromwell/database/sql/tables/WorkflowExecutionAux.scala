package cromwell.database.sql.tables

import java.sql.Clob

@deprecated("Olde Worlde Databasee Tablee", "0.21")
case class WorkflowExecutionAux
(
  workflowExecutionId: Int,
  wdlSource: Clob,
  jsonInputs: Clob,
  workflowOptions: Clob,
  workflowExecutionAuxId: Option[Int] = None
)
