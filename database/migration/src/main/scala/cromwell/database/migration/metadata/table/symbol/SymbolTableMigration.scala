package cromwell.database.migration.metadata.table.symbol

import java.sql.{PreparedStatement, ResultSet}

import cromwell.core.simpleton.WdlValueSimpleton._
import cromwell.database.migration.WdlTransformation
import cromwell.database.migration.custom.BatchedTaskChange
import wdl4s.WdlExpression
import wdl4s.types.WdlType
import wdl4s.values.WdlValue

import scala.util.{Failure, Success, Try}

object SymbolTableMigration {
  val NbRowsQuery =
    """
      |SELECT MAX(TMP_SYMBOL_ID) AS symbol_count
      |FROM TMP_SYMBOL;
    """.stripMargin
}

trait SymbolTableMigration extends BatchedTaskChange {
  import cromwell.database.migration.WdlTransformation._

  override val readCountQuery = SymbolTableMigration.NbRowsQuery

  override val readBatchQuery = """
      |SELECT
      |    WORKFLOW_EXECUTION_UUID,
      |    SYMBOL_NAME,
      |    SYMBOL_SCOPE,
      |    SYMBOL_INDEX,
      |    SYMBOL_ATTEMPT,
      |    WDL_TYPE,
      |    WDL_VALUE
      |   FROM TMP_SYMBOL
      |   WHERE TMP_SYMBOL_ID >= ? AND TMP_SYMBOL_ID < ?;
    """.stripMargin

  override val migrateBatchQuery = MetadataStatement.InsertSql

  /**
    * Migrate a row to the metadata table
    */
  override def migrateBatchRow(row: ResultSet, statement: PreparedStatement): Int = {
    // Try to coerce the value to a WdlValue
    val value = for {
      wdlType <- Try(WdlType.fromWdlString(row.getString("WDL_TYPE")))
      inflated <- row.getString("WDL_VALUE") match {
        case null => Success("") // Empty Strings are null in the DB
        case nonNull => inflate(row.getString("WDL_VALUE"))
      }
    } yield WdlTransformation.coerceStringToWdl(inflated, wdlType)

    val workflowUuid = row.getString("WORKFLOW_EXECUTION_UUID")
    val symbolName = row.getString("SYMBOL_NAME")
    val symbolScope = row.getString("SYMBOL_SCOPE")
    // getInt returns 0 if value is null so use getString instead and convert after
    val symbolIndex = Option(row.getString("SYMBOL_INDEX")) map { _.toInt }
    val symbolAttempt = Option(row.getString("SYMBOL_ATTEMPT")) map { _.toInt }

    value match {
      case Success(wdlValue) =>
        processSymbol(statement, workflowUuid, symbolName, symbolScope, symbolIndex, symbolAttempt, wdlValue)
      case Failure(f) =>
        logger.error(
          s"""Could not parse symbol of type ${row.getString("WDL_TYPE")}
              |for Workflow $workflowUuid - Call $symbolScope:$symbolIndex""".stripMargin, f)
        0
    }
  }

  def processSymbol(statement: PreparedStatement,
                    workflowUuid: String,
                    symbolName: String,
                    symbolScope: String,
                    symbolIndex: Option[Int],
                    symbolAttempt: Option[Int],
                    wdlValue: WdlValue): Int

  /**
    * Add all necessary statements to the batch for the provided WdlValue.
    */
  protected final def addWdlValue(metadataKey: String, wdlValue: WdlValue, metadataStatementForCall: MetadataStatement): Int = {
    wdlValue match {
        // simplify doesn't handle WdlExpression
      case expr: WdlExpression =>
        metadataStatementForCall.addKeyValue(metadataKey, expr.valueString)
        1
      case value =>
        val simplified = value.simplify(metadataKey)
        simplified foreach { simpleton =>
          metadataStatementForCall.addKeyValue(simpleton.simpletonKey, simpleton.simpletonValue)
        }
        simplified.size
    }
  }
}
