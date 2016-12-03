package cromwell.database.migration.metadata.table.symbol

import java.sql.PreparedStatement

import wdl4s.values._

class InputSymbolTableMigration extends SymbolTableMigration {

  override def processSymbol(statement: PreparedStatement,
                             workflowUuid: String,
                             symbolName: String,
                             symbolScope: String,
                             symbolIndex: Option[Int],
                             symbolAttempt: Option[Int],
                             wdlValue: WdlValue): Int = {

    (symbolIndex, symbolAttempt) match {
      case (Some(index) , Some(attempt)) =>
        // Call scoped
        val metadataStatementForCall = new MetadataStatementForCall(statement,
          workflowUuid,
          symbolScope,
          index,
          attempt
        )

        addWdlValue(s"inputs:$symbolName", wdlValue, metadataStatementForCall)
      case (None, None) if !symbolScope.contains('.') =>
        val metadataStatementForWorkflow = new MetadataStatementForWorkflow(statement, workflowUuid)
        addWdlValue(s"inputs:$symbolScope.$symbolName", wdlValue, metadataStatementForWorkflow)
      case _ =>
        0
    }
  }

  override def getConfirmationMessage: String = "Inputs from Symbol Table migration complete."
}
