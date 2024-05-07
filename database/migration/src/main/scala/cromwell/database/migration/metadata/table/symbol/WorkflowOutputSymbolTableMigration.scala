package cromwell.database.migration.metadata.table.symbol

import java.sql.PreparedStatement

import wom.values.WomValue

class WorkflowOutputSymbolTableMigration extends SymbolTableMigration {

  override def processSymbol(statement: PreparedStatement,
                             workflowUuid: String,
                             symbolName: String,
                             symbolScope: String,
                             symbolIndex: Option[Int],
                             symbolAttempt: Option[Int],
                             womValue: WomValue
  ): Int = {
    val metadataStatementForWorkflow = new MetadataStatementForWorkflow(statement, workflowUuid)
    addWdlValue(s"outputs:$symbolScope.$symbolName", womValue, metadataStatementForWorkflow)
  }

  override def getConfirmationMessage: String = "Workflow outputs from Symbol Table migration complete."
}
