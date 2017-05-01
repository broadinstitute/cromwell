package cromwell.database.migration.failuremetadata

import java.sql.{PreparedStatement, ResultSet}
import cromwell.database.migration.custom.BatchedTaskChange

/**
  * Turns individual '| failures[0] | "blah" |' entries into '| failures[0]:message | "blah" |' and '| failures[0]:causedBy[] | NULL |'
  */
class ExpandSingleFailureStrings extends BatchedTaskChange {
  val tableName = "METADATA_ENTRY"
  val primaryKeyColumn = "METADATA_JOURNAL_ID"
  val workflowIdColumn = "WORKFLOW_EXECUTION_UUID"
  val metadataKeyColumn = "METADATA_KEY"
  val callFqnColumn = "CALL_FQN"
  val jobScatterIndexColumn = "JOB_SCATTER_INDEX"
  val retryAttemptColumn = "JOB_RETRY_ATTEMPT"
  val metadataTimestampColumn = "METADATA_TIMESTAMP"

  val fixableFailureMessageFilter = "METADATA_KEY REGEXP '^failures\\\\[[0-9]*\\\\]$'"

  override def readCountQuery: String = s"SELECT MAX($primaryKeyColumn) FROM $tableName;"

  /**
    * Find metadata keys in this range which can be fixed up.
    */
  override val readBatchQuery: String =
    s"""|SELECT $primaryKeyColumn, $workflowIdColumn, $metadataKeyColumn, $callFqnColumn, $jobScatterIndexColumn, $retryAttemptColumn, $metadataTimestampColumn
        | FROM $tableName
        | WHERE $primaryKeyColumn >= ? AND $primaryKeyColumn < ?
        | AND $fixableFailureMessageFilter
        |""".stripMargin

  override val migrateBatchQueries = List(
    // Modify the existing element, changing failures[x] to failures[x]:message
    s"UPDATE $tableName SET $metadataKeyColumn = ? WHERE $primaryKeyColumn = ?",
    // Insert a new element, to contain the empty causedBy[]
    s"""INSERT INTO $tableName ($workflowIdColumn, $metadataKeyColumn, $callFqnColumn, $jobScatterIndexColumn, $retryAttemptColumn, $metadataTimestampColumn)
       | VALUES (?, ?, ?, ?, ?, ?)
     """.stripMargin
  )

  override def migrateBatchRow(readRow: ResultSet, migrateStatements: List[PreparedStatement]): Int = {
    val (modifyStatement :: insertStatement :: Nil) = migrateStatements

    val rowId = readRow.getInt(primaryKeyColumn)
    val workflowId = readRow.getString(workflowIdColumn)
    val oldKey = readRow.getString(metadataKeyColumn)
    val callFqn = readRow.getString(callFqnColumn)
    val jobScatterIndex = readRow.getObject(jobScatterIndexColumn)
    val retryAttempt = readRow.getObject(retryAttemptColumn)
    val metadataTimestamp = readRow.getObject(metadataTimestampColumn)

    val modifiedKey = oldKey + ":message"
    val insertKey = oldKey + ":causedBy[]"

    modifyStatement.setString(1, modifiedKey)
    modifyStatement.setInt(2, rowId)
    modifyStatement.addBatch()

    insertStatement.setString(1, workflowId)
    insertStatement.setString(2, insertKey)
    insertStatement.setString(3, callFqn)
    insertStatement.setObject(4, jobScatterIndex)
    insertStatement.setObject(5, retryAttempt)
    insertStatement.setObject(6, metadataTimestamp)
    insertStatement.addBatch()

    2
  }
}
