package cromwell.database.migration.workflowoptions

import cromwell.core.WorkflowOptions

/**
  * Clear the values from encrypted keys in METADATA_ENTRY.
  */
class ClearMetadataEntryWorkflowOptions extends WorkflowOptionsChange {
  override val tableName = "METADATA_ENTRY"
  override val primaryKeyColumn = "METADATA_JOURNAL_ID"
  override val workflowOptionsColumn = "METADATA_VALUE"
  override val additionalReadBatchFilters = "AND METADATA_KEY = 'submittedFiles:options'"

  override def migrateWorkflowOptions(workflowOptions: WorkflowOptions) = workflowOptions.clearEncryptedValues
}
