package cromwell.database.migration.restart.table

import cromwell.database.migration.workflowoptions.WorkflowOptionsRenaming._
import cromwell.database.migration.restart.table.RenameWorkflowOptionKeysMigration._
import liquibase.database.jvm.JdbcConnection
import spray.json._

object RenameWorkflowOptionKeysMigration {
  private val QueryWorkflowStore = " SELECT WORKFLOW_STORE_ID, WORKFLOW_OPTIONS FROM WORKFLOW_STORE "

  private val UpdateWorkflowStore = " UPDATE WORKFLOW_STORE SET WORKFLOW_OPTIONS = ? WHERE WORKFLOW_STORE_ID = ? "
}


class RenameWorkflowOptionKeysMigration extends AbstractRestartMigration {

  override protected def description: String = "Workflow option renaming"

  override protected def doMigration(connection: JdbcConnection): Unit = {
    val query = connection.createStatement()
    lazy val insert = connection.prepareStatement(UpdateWorkflowStore)
    query.execute(QueryWorkflowStore)
    val results = query.getResultSet
    while (results.next()) {
      val options = results.getString(2)
      if (!results.wasNull()) {
        val optionsJson = options.parseJson
        val newOptionsJson = optionsJson match {
          case JsObject(fields) => JsObject(fields map renameOptionKeys)
          case other => other // There really shouldn't be workflow options of other types, but if there are pass them through.
        }

        insert.setString(1, newOptionsJson.prettyPrint)
        insert.setInt(2, results.getInt(1))
        insert.execute()
      }
    }
  }
}
