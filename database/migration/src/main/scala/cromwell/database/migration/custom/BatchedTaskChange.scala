package cromwell.database.migration.custom

import java.sql.{PreparedStatement, ResultSet}

import liquibase.database.jvm.JdbcConnection
import liquibase.exception.CustomChangeException

/**
  * Runs a migration as a series of batches.
  */
trait BatchedTaskChange extends MigrationTaskChange {
  /**
    * Returns sql to retrieve the maximum primary key for the table.
    *
    * Example:
    * {{{
    *   SELECT MAX([PRIMARY_KEY])
    *   FROM [TABLE];
    * }}}
    */
  def readCountQuery: String

  /**
    * Returns sql to retrieve rows to be passed to migrateBatchRow, batching on a primary key between the half-open
    * primary key range [start, stop).
    *
    * Example:
    * {{{
    *   SELECT [COLUMNS]
    *   FROM [TABLE]
    *   WHERE [PRIMARY_KEY] >= ? AND [PRIMARY_KEY] < ?;
    * }}}
    */
  def readBatchQuery: String

  /**
    * Used to prepare the statement that will be passed repeatedly into migrateBatchRow.
    *
    * Example:
    * {{{
    *   UPDATE [TABLE]
    *   SET [COLUMNS]
    *   WHERE [PRIMARY_KEY] = ?;
    * }}}
    *
    * Example:
    * {{{
    *   INSERT INTO [TABLE]
    *   SET [COLUMNS];
    * }}}
    */
  def migrateBatchQuery: String

  /**
    * Migrate a row.
    *
    * Read the values from readRow, update the values, set the updated values on the migrateStatement, and then call
    * migrateStatement.addBatch(). Return the (estimated) number of rows to be written by this batch.
    *
    * @param readRow          The row to migrate
    * @param migrateStatement The statement to add a new migrated row
    * @return The number of rows updated
    */
  def migrateBatchRow(readRow: ResultSet, migrateStatement: PreparedStatement): Int

  /**
    * Specify the size of a "page".
    * For databases with a very large number of rows, selecting all the rows at once can generate a variety of problems.
    * In order to avoid any issue, the selection is paginated. This value sets how many rows should be retrieved and
    * processed at a time, before asking for the next page.
    */
  private val readBatchSize = config.getInt("database.migration.read-batch-size")

  /**
    * To keep the size of the insert batch from growing out of control we monitor its size and execute/commit when it
    * reaches or exceeds writeBatchSize.
    */
  private val writeBatchSize = config.getInt("database.migration.write-batch-size")

  override def migrate(connection: JdbcConnection) = {

    logger.info(s"Running migration $migrationName with a read batch size of " +
      s"$readBatchSize and a write batch size of $writeBatchSize")

    /*
      * Keep count of the size of the batch.
      *
      * @see writeBatchSize
      */
    var batchMigrationCounter: Int = 0

    val readCount = getReadCount(connection)

    // So we can display progress
    val pageCount = Math.max(readCount / readBatchSize, 1)

    val readBatchStatement = connection.prepareStatement(readBatchQuery)
    val migrateBatchStatement = connection.prepareStatement(migrateBatchQuery)

    val paginator = new QueryPaginator(readBatchStatement, readBatchSize, readCount)

    // Loop over pages
    paginator.zipWithIndex foreach {
      case (resultBatch, page) =>
        // Loop over rows in page
        new ResultSetIterator(resultBatch).zipWithIndex foreach {
          case (row, idx) =>
            batchMigrationCounter += migrateBatchRow(row, migrateBatchStatement)
            // batchMigrationCounter can actually be bigger than writeBatchSize as wdlValues are processed atomically,
            // so this is a best effort
            if (batchMigrationCounter >= writeBatchSize) {
              migrateBatchStatement.executeBatch()
              connection.commit()
              batchMigrationCounter = 0
            }
        }

        resultBatch.close()

        val progress = Math.min((page + 1) * 100 / pageCount, 100)
        logger.info(s"[$migrationName] $progress%")
    }

    if (batchMigrationCounter != 0) {
      migrateBatchStatement.executeBatch()
      connection.commit()
    }
  }

  private def getReadCount(connection: JdbcConnection): Int = {
    val readCountResultSet = connection.createStatement().executeQuery(readCountQuery)

    if (readCountResultSet.next()) {
      readCountResultSet.getInt(1)
    } else {
      throw new CustomChangeException(s"Could not find max value for pagination from sql:\n$readCountQuery")
    }
  }
}
