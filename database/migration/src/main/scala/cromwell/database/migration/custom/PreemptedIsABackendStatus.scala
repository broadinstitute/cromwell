package cromwell.database.migration.custom

import java.sql.{PreparedStatement, ResultSet}

/**
  * Makes sure that it's executionStatuses that are 'Failed' and backendStatuses can be 'Preempted'
  */
class PreemptedIsABackendStatus extends BatchedTaskChange {
  val tableName = "METADATA_ENTRY"
  val primaryKeyColumn = "METADATA_JOURNAL_ID"
  val workflowIdColumn = "WORKFLOW_EXECUTION_UUID"
  val metadataKeyColumn = "METADATA_KEY"
  val callFqnColumn = "CALL_FQN"
  val jobScatterIndexColumn = "JOB_SCATTER_INDEX"
  val retryAttemptColumn = "JOB_RETRY_ATTEMPT"
  val timestampColumn = "METADATA_TIMESTAMP"
  val metadataValueColumn = "METADATA_VALUE"
  val metadataValueTypeColumn = "METADATA_VALUE_TYPE"

  val messageFilter = "METADATA_KEY = 'executionStatus' && METADATA_VALUE = 'Preempted'"

  override def readCountQuery = s"SELECT MAX($primaryKeyColumn) FROM $tableName"

  /**
    * Find metadata keys in this range which are duplicates of earlier keys (for the same workflow).
    */
  override def readBatchQuery =
    s"""|SELECT $primaryKeyColumn, $workflowIdColumn, $callFqnColumn, $jobScatterIndexColumn, $retryAttemptColumn, $timestampColumn
        | FROM $tableName AS t1
        | WHERE t1.$primaryKeyColumn >= ? AND t1.$primaryKeyColumn < ?
        | AND $messageFilter
        |""".stripMargin

  override def migrateBatchQueries = List(
    s"UPDATE $tableName SET $metadataKeyColumn = 'backendStatus' WHERE $primaryKeyColumn = ?",
    s"INSERT INTO $tableName ($workflowIdColumn, $callFqnColumn, $jobScatterIndexColumn, $retryAttemptColumn, $timestampColumn, $metadataKeyColumn, $metadataValueColumn, $metadataValueTypeColumn) " +
      s"VALUES (?, ?, ?, ?, ?, 'executionStatus', 'RetryableFailure', 'string')"
  )

  override def migrateBatchRow(readRow: ResultSet, migrateStatements: List[PreparedStatement]): Int = {
    val rowId = readRow.getInt(primaryKeyColumn)
    val workflowId = readRow.getString(workflowIdColumn)
    val callFqn = readRow.getString(callFqnColumn)
    val jobScatterIndex = readRow.getObject(jobScatterIndexColumn)
    val retryAttempt = readRow.getObject(retryAttemptColumn)
    val metadataTimestamp = readRow.getObject(timestampColumn)

    // First, update the executionStatus to be a backendStatus
    val updateStatement = migrateStatements.head
    updateStatement.setInt(1, rowId)
    updateStatement.addBatch()

    // Second, insert the new 'RetryableFailure' executionStatus
    val insertStatement = migrateStatements(1)
    insertStatement.setString(1, workflowId)
    insertStatement.setString(2, callFqn)
    insertStatement.setObject(3, jobScatterIndex)
    insertStatement.setObject(4, retryAttempt)
    insertStatement.setObject(5, metadataTimestamp)
    insertStatement.addBatch()

    // We added two queries here, so that's what we return:
    2
  }
}
