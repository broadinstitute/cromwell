package cromwell.database.migration.metadata.table.symbol

import java.sql.PreparedStatement

import wdl4s.values._

class WorkflowOutputSymbolTableMigration extends SymbolTableMigration {

  override def processSymbol(statement: PreparedStatement,
                             workflowUuid: String,
                             symbolName: String,
                             symbolScope: String,
                             symbolIndex: Option[Int],
                             symbolAttempt: Option[Int],
                             wdlValue: WdlValue): Int = {
    val metadataStatementForWorkflow = new MetadataStatementForWorkflow(statement, workflowUuid)
    addWdlValue(s"outputs:$symbolScope.$symbolName", wdlValue, metadataStatementForWorkflow)
  }

  override def getConfirmationMessage: String = "Workflow outputs from Symbol Table migration complete."
}
