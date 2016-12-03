package cromwell.database.migration.workflowoptions

import cromwell.core.WorkflowOptions

/**
  * Encrypt the values for encrypted keys in WORKFLOW_STORE_ENTRY.
  */
class EncryptWorkflowStoreEntryWorkflowOptions extends WorkflowOptionsChange {
  override val tableName = "WORKFLOW_STORE_ENTRY"
  override val primaryKeyColumn = "WORKFLOW_STORE_ENTRY_ID"
  override val workflowOptionsColumn = "WORKFLOW_OPTIONS"

  override def migrateWorkflowOptions(workflowOptions: WorkflowOptions) = workflowOptions.asPrettyJson
}
