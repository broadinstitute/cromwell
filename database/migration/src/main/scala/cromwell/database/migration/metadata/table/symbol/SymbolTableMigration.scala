package cromwell.database.migration.metadata.table.symbol

import java.sql.{PreparedStatement, ResultSet}

import com.typesafe.config.ConfigFactory
import cromwell.core.simpleton.WdlValueSimpleton._
import cromwell.database.migration.WdlTransformation
import liquibase.change.custom.CustomTaskChange
import liquibase.database.Database
import liquibase.database.jvm.JdbcConnection
import liquibase.exception.{CustomChangeException, ValidationErrors}
import liquibase.resource.ResourceAccessor
import org.slf4j.LoggerFactory
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

trait SymbolTableMigration extends CustomTaskChange {
  import SymbolTableMigration._
  import cromwell.database.migration.WdlTransformation._

  // Nb of rows to retrieve / process in a batch
  val config = ConfigFactory.load

  /**
   * Specify the size of a "page".
   * For databases with a very large number of symbols, selecting all the rows at once can generate a variety of problems.
   * In order to avoid any issue, the selection is paginated. This value sets how many rows should be retrieved and processed at a time, before asking for the next chunk.
   */
  val readBatchSize = config.getInt("database.migration.read-batch-size")

  /**
   * Because a symbol row can contain any arbitrary wdl value, the amount of metadata rows to insert from a single symbol row can vary from 1 to several thousands (or more).
   * To keep the size of the insert batch from growing out of control we monitor its size and execute/commit when it reaches or exceeds writeBatchSize.
   */
  val writeBatchSize = config.getInt("database.migration.write-batch-size")

  val logger = LoggerFactory.getLogger("LiquibaseMetadataMigration")

  override def execute(database: Database): Unit = {
    try {
      val dbConn = database.getConnection.asInstanceOf[JdbcConnection]
      val autoCommit = dbConn.getAutoCommit
      dbConn.setAutoCommit(false)
      migrate(dbConn)
      dbConn.setAutoCommit(autoCommit)
    } catch {
      case t: CustomChangeException => throw t
      case t: Throwable => throw new CustomChangeException(s"Could not apply migration script for metadata at ${getClass.getSimpleName}", t)
    }
  }

  def tmpSymbolPaginatedStatement(connection: JdbcConnection): PreparedStatement = connection.prepareStatement("""
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
    """.stripMargin)

  private def migrate(connection: JdbcConnection) = {
    logger.info(s"Running migration with a read batch size of $readBatchSize and a write batch size of $writeBatchSize")

    /**
     * Keep count of the size of the batch.
      *
      * @see writeBatchSize
     */
    var insertsCounter: Int = 0

    // Find the max row id in the TMP_SYMBOL table
    val tmpSymbolCountRS = connection.createStatement().executeQuery(NbRowsQuery)

    if (tmpSymbolCountRS.next()) {
      val tmpSymbolCount = tmpSymbolCountRS.getInt("symbol_count")

      // So we can display progress
      val nbPages = Math.max(tmpSymbolCount / readBatchSize, 1)

      val paginator = new QueryPaginator(tmpSymbolPaginatedStatement(connection), readBatchSize, tmpSymbolCount)
      val metadataInsertStatement = MetadataStatement.makeStatement(connection)

      // Loop over pages
      paginator.zipWithIndex foreach {
        case (resultBatch, page) =>
          // Loop over rows in page
          new ResultSetIterator(resultBatch).zipWithIndex foreach {
            case (row, idx) =>
              insertsCounter += migrateRow(connection, metadataInsertStatement, row, idx)
              // insertsCounter can actually be bigger than writeBatchSize as wdlValues are processed atomically, so this is a best effort
              if (insertsCounter >= writeBatchSize) {
                metadataInsertStatement.executeBatch()
                connection.commit()
                insertsCounter = 0
              }
          }

          resultBatch.close()

          val progress = Math.min((page + 1) * 100 / nbPages, 100)
          logger.info(s"[${getClass.getSimpleName}] $progress%")
      }

      if (insertsCounter != 0) {
        metadataInsertStatement.executeBatch()
        connection.commit()
      }
    } else {
      throw new CustomChangeException("Could not find max value of symbol id for pagination")
    }
  }

  /**
    * Migrate a row to the metadata table
    */
  protected def migrateRow(connection: JdbcConnection, statement: PreparedStatement, row: ResultSet, idx: Int): Int = {
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
        processSymbol(statement, idx, workflowUuid, symbolName, symbolScope, symbolIndex, symbolAttempt, wdlValue)
      case Failure(f) =>
        logger.error(
          s"""Could not parse symbol of type ${row.getString("WDL_TYPE")}
              |for Workflow $workflowUuid - Call $symbolScope:$symbolIndex""".stripMargin, f)
        0
    }
  }

  def processSymbol(statement: PreparedStatement,
                    idx: Int,
                    workflowUuid: String,
                    symbolName: String,
                    symbolScope: String,
                    symbolIndex: Option[Int],
                    symbolAttempt: Option[Int],
                    wdlValue: WdlValue): Int

  override def setUp(): Unit = ()

  override def validate(database: Database): ValidationErrors = new ValidationErrors

  override def setFileOpener(resourceAccessor: ResourceAccessor): Unit = {}

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
