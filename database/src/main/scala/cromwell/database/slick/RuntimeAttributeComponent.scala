package cromwell.database.slick

import cromwell.database.obj.RuntimeAttribute

trait RuntimeAttributeComponent {
  this: DriverComponent with ExecutionComponent with WorkflowExecutionComponent =>

  import driver.api._

  class RuntimeAttributes(tag: Tag) extends Table[RuntimeAttribute](tag, "RUNTIME_ATTRIBUTES") {
    def runtimeAttributeId = column[Int]("RUNTIME_ATTRIBUTE_ID", O.PrimaryKey, O.AutoInc)
    def executionId = column[Int]("EXECUTION_ID")
    def name = column[String]("ATTRIBUTE_NAME")
    def value = column[String]("ATTRIBUTE_VALUE")
    override def * = (executionId, name, value, runtimeAttributeId.?) <> (RuntimeAttribute.tupled, RuntimeAttribute.unapply)
    def execution = foreignKey("FK_RUNTIME_ATTRIBUTE_EXECUTION", executionId, executions)(_.executionId)
    def uniqueKey = index("UK_RUNTIME_ATTRIBUTE", (executionId, name), unique = true)
  }

  protected val runtimeAttributes = TableQuery[RuntimeAttributes]
  val runtimeAttributeIdsAutoInc = runtimeAttributes returning runtimeAttributes.map(_.runtimeAttributeId)

  val runtimeAttributeValuesByExecutionIdAndName = Compiled(
    (executionId: Rep[Int], name: Rep[String]) => for {
      runtimeAttribute <- runtimeAttributes
      if runtimeAttribute.name === name
      if runtimeAttribute.executionId === executionId
    } yield runtimeAttribute.value)

  val runtimeAttributesByWorkflowUUID = Compiled(
    (workflowUUID: Rep[String]) => for {
      workflow <- workflowExecutions
      if workflow.workflowExecutionUuid === workflowUUID
      execution <- executions
      if execution.workflowExecutionId === workflow.workflowExecutionId
      runtimeAttribute <- runtimeAttributes
      if runtimeAttribute.executionId === execution.executionId
    } yield (execution.callFqn, execution.index, execution.attempt, runtimeAttribute.name, runtimeAttribute.value))
}
