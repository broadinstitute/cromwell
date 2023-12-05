package cromwell.database.migration.metadata.table.symbol

import java.sql.PreparedStatement

import wom.values.WomValue

class CallOutputSymbolTableMigration extends SymbolTableMigration {
  override def processSymbol(statement: PreparedStatement,
                             workflowUuid: String,
                             symbolName: String,
                             symbolScope: String,
                             symbolIndex: Option[Int],
                             symbolAttempt: Option[Int],
                             womValue: WomValue
  ): Int =
    (symbolIndex, symbolAttempt) match {
      case (Some(index), Some(attempt)) =>
        val metadataStatementForCall =
          new MetadataStatementForCall(statement, workflowUuid, symbolScope, index, attempt)

        addWdlValue(s"outputs:$symbolName", womValue, metadataStatementForCall)
      case _ =>
        logger.warn(
          s"Found output without index or attempt: [$workflowUuid] $symbolScope - $symbolName:$symbolIndex:$symbolAttempt"
        )
        0
    }

  override def getConfirmationMessage: String = "Call outputs from Symbol Table migration complete."
}
