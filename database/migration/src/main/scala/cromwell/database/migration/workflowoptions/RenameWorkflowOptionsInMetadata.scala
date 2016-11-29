package cromwell.database.migration.workflowoptions

import java.sql.{PreparedStatement, ResultSet}

import cromwell.database.migration.custom.BatchedTaskChange
import cromwell.database.migration.workflowoptions.WorkflowOptionsRenaming._
import spray.json.{JsObject, _}

class RenameWorkflowOptionsInMetadata extends BatchedTaskChange {
  val tableName = "METADATA_ENTRY"
  val primaryKeyColumn = "METADATA_JOURNAL_ID"
  val workflowOptionsColumn = "METADATA_VALUE"
  val additionalReadBatchFilters = "AND METADATA_KEY = 'submittedFiles:options'"

  override def readCountQuery = s"SELECT MAX($primaryKeyColumn) FROM $tableName;"

  override def readBatchQuery =
    s"""|SELECT $primaryKeyColumn, $workflowOptionsColumn
        |  FROM $tableName
        |  WHERE $primaryKeyColumn >= ? AND $primaryKeyColumn < ? $additionalReadBatchFilters;
        |""".stripMargin

  override def migrateBatchQuery = s"UPDATE $tableName SET $workflowOptionsColumn = ? WHERE $primaryKeyColumn = ?;"

  override def migrateBatchRow(readRow: ResultSet, migrateStatement: PreparedStatement): Int = {
    val rowId = readRow.getInt(1)
    
    val migratedJson = readRow.getString(2).parseJson match {
      case JsObject(fields) => JsObject(fields map renameOptionKeys)
      case other => other
    }

    migrateStatement.setString(1, migratedJson.prettyPrint)
    migrateStatement.setInt(2, rowId)
    migrateStatement.addBatch()
    1
  }
}
