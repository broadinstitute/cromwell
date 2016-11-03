package cromwell.database.migration.workflowoptions

import java.sql.{PreparedStatement, ResultSet}

import cromwell.core.WorkflowOptions
import cromwell.database.migration.custom.BatchedTaskChange
import liquibase.database.Database

import scala.util.{Failure, Success}

/**
  * Edits the workflow options stored in a table.
  */
trait WorkflowOptionsChange extends BatchedTaskChange {
  /** @return name of the table */
  def tableName: String

  /** @return primary key of the table */
  def primaryKeyColumn: String

  /** @return column storing the workflow options */
  def workflowOptionsColumn: String

  /** @return any additional filters to add to the where clause, starting with "AND ..." */
  def additionalReadBatchFilters: String = ""

  /**
    * Takes in the workflow options and returns the edited version as a json string.
    *
    * @param workflowOptions workflow options object
    * @return edited workflow object json
    */
  def migrateWorkflowOptions(workflowOptions: WorkflowOptions): String

  override def execute(database: Database): Unit = {
    val configPath = "workflow-options.encrypted-fields"
    if (config.hasPath(configPath) && !config.getStringList(configPath).isEmpty) {
      super.execute(database)
    }
  }

  override def readCountQuery = s"SELECT MAX($primaryKeyColumn) FROM $tableName;"

  override def readBatchQuery =
    s"""|SELECT $primaryKeyColumn, $workflowOptionsColumn
        |  FROM $tableName
        |  WHERE $primaryKeyColumn >= ? AND $primaryKeyColumn < ? $additionalReadBatchFilters;
        |""".stripMargin

  override def migrateBatchQuery = s"UPDATE $tableName SET $workflowOptionsColumn = ? WHERE $primaryKeyColumn = ?;"

  override def migrateBatchRow(readRow: ResultSet, migrateStatement: PreparedStatement): Int = {
    val rowId = readRow.getInt(1)
    val workflowOptionsJson = readRow.getString(2)
    WorkflowOptions.fromJsonString(workflowOptionsJson) match {
      case Success(workflowOptions) =>
        val migratedJson = migrateWorkflowOptions(workflowOptions)
        migrateStatement.setString(1, migratedJson)
        migrateStatement.setInt(2, rowId)
        migrateStatement.addBatch()
        1
      case Failure(exception) =>
        logger.error(
          s"Unable to process $tableName pk $rowId\njson:\n$workflowOptionsJson", exception)
        0
    }
  }

}
