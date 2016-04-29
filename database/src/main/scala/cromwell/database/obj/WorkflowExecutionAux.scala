package cromwell.database.obj

import java.sql.Clob

case class WorkflowExecutionAux
(
  workflowExecutionId: Int,
  wdlSource: Clob,
  jsonInputs: Clob,
  workflowOptions: Clob,
  workflowExecutionAuxId: Option[Int] = None
)
