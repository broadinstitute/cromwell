package cromwell.engine.db.slick

// TODO switch to UUID on WorkflowExecution rather than synthetic PK, update this to reference that.
// TODO probably don't need a synthetic PK on this table as the WorkflowExecution ID serves that purpose.
case class WorkflowExecutionAux
(
  workflowExecutionId: Int,
  wdlSource: String,
  jsonInputs: String,
  workflowExecutionAuxId: Option[Int] = None
)

trait WorkflowExecutionAuxComponent {
  this: DriverComponent with WorkflowExecutionComponent =>

  import driver.api._

  class WorkflowExecutionAuxes(tag: Tag) extends Table[WorkflowExecutionAux](tag, "WORKFLOW_EXECUTION_AUX") {
    def workflowExecutionAuxId = column[Int]("WORKFLOW_EXECUTION_AUX_ID", O.PrimaryKey, O.AutoInc)
    def workflowExecutionId = column[Int]("WORKFLOW_EXECUTION_ID")
    def wdlSource = column[String]("WDL_SOURCE")
    def jsonInputs = column[String]("JSON_INPUTS")

    override def * = (workflowExecutionId, wdlSource, jsonInputs, workflowExecutionAuxId.?) <>
      (WorkflowExecutionAux.tupled, WorkflowExecutionAux.unapply)

    def workflow = foreignKey(
      "FK_WE_AUX_WORKFLOW_EXECUTION_ID", workflowExecutionId, workflowExecutions)(_.workflowExecutionId)

  }

  val workflowExecutionAuxes = TableQuery[WorkflowExecutionAuxes]

  val workflowExecutionAuxesAutoInc = workflowExecutionAuxes returning workflowExecutionAuxes.
    map(_.workflowExecutionAuxId) into ((a, id) => a.copy(workflowExecutionAuxId = Some(id)))

}
