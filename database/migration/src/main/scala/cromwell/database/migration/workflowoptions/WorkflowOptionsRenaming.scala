package cromwell.database.migration.workflowoptions

import spray.json._

object WorkflowOptionsRenaming {

  private val RenamedOptionKeys = Map(
    "defaultRuntimeOptions" -> "default_runtime_attributes",
    "workflowFailureMode" -> "workflow_failure_mode",
    "workflow_log_dir" -> "final_workflow_log_dir",
    "outputs_path" -> "final_workflow_outputs_dir",
    "call_logs_dir" -> "final_call_logs_dir"
  )

  def renameOptionKeys(field: JsField): JsField = {
    field match {
      case (oldName, value) if RenamedOptionKeys.contains(oldName) => RenamedOptionKeys(oldName) -> value
      case noop => noop
    }
  }
}
