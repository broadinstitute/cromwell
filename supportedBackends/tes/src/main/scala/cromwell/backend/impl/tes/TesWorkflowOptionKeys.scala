package cromwell.backend.impl.tes

object TesWorkflowOptionKeys {
  // Communicates to the TES server which identity the tasks should execute as
  val WorkflowExecutionIdentity = "workflow_execution_identity"
  val DataAccessIdentity = "data_access_identity"
}
