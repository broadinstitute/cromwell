package cromwell.database.migration.failuremetadata

import java.sql.{PreparedStatement, ResultSet}

import cromwell.database.migration.custom.BatchedTaskChange

import scala.util.Random

/**
  * Turns 'causedBy' objects into lists in the metadata.
  */
class DeduplicateFailureMessageIds extends BatchedTaskChange {
  val tableName = "METADATA_ENTRY"
  val primaryKeyColumn = "METADATA_JOURNAL_ID"
  val workflowIdColumn = "WORKFLOW_EXECUTION_UUID"
  val metadataKeyColumn = "METADATA_KEY"
  val callFqnColumn = "CALL_FQN"
  val jobScatterIndexColumn = "JOB_SCATTER_INDEX"
  val retryAttemptColumn = "JOB_RETRY_ATTEMPT"
  val contentEqualityCheck = List(workflowIdColumn, metadataKeyColumn, callFqnColumn, jobScatterIndexColumn, retryAttemptColumn)
    .map(s => s"(t2.$s = t1.$s OR (t2.$s IS NULL AND t1.$s IS NULL))").mkString(" AND ")

  val fixableFailureMessageFilter = "METADATA_KEY LIKE '%failures[%]%:message'"

  override def readCountQuery = s"SELECT MAX($primaryKeyColumn) FROM $tableName;"

  /**
    * Find metadata keys in this range which are duplicates of earlier keys (for the same workflow).
    */
  override def readBatchQuery =
    s"""|SELECT $primaryKeyColumn, $workflowIdColumn, $metadataKeyColumn
        | FROM $tableName AS t1
        | WHERE t1.$primaryKeyColumn >= ? AND t1.$primaryKeyColumn < ?
        | AND $fixableFailureMessageFilter
        | AND EXISTS (SELECT *
        |              FROM $tableName AS t2
        |              WHERE $contentEqualityCheck
        |                AND t2.$primaryKeyColumn < t1.$primaryKeyColumn
        |            )
        |""".stripMargin

  override def migrateBatchQuery = s"UPDATE $tableName SET $metadataKeyColumn = ? WHERE $primaryKeyColumn = ?;"

  override def migrateBatchRow(readRow: ResultSet, migrateStatement: PreparedStatement): Int = {
    val rowId = readRow.getInt(1)
    val oldKey = readRow.getString(3)
    val newRandomInt = Random.nextInt(Int.MaxValue)
    val newKey = oldKey.replaceFirst("failures\\[\\d*\\]", s"failures[$newRandomInt]")

    migrateStatement.setString(1, newKey)
    migrateStatement.setInt(2, rowId)
    migrateStatement.addBatch()
    1
  }
}
