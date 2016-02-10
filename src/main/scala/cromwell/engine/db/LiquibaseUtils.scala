package cromwell.engine.db

import java.sql.Connection

import liquibase.database.jvm.{HsqlConnection, JdbcConnection}
import liquibase.database.{Database, DatabaseConnection, DatabaseFactory}
import liquibase.diff.compare.CompareControl
import liquibase.diff.{DiffGeneratorFactory, DiffResult}
import liquibase.resource.ClassLoaderResourceAccessor
import liquibase.{Contexts, LabelExpression, Liquibase}
import org.hsqldb.persist.HsqlDatabaseProperties

import scala.collection.JavaConverters._

object LiquibaseUtils {
  val ChangeLogDir = "src/main/migrations/"
  val DefaultChangeLog = "changelog.xml"
  val DefaultContexts = new Contexts()
  val DefaultLabelExpression = new LabelExpression()

  /**
    * Updates a liquibase schema to the latest version.
    *
    * @param jdbcConnection A jdbc connection to the database.
    */
  def updateSchema(jdbcConnection: Connection): Unit = {
    val liquibaseConnection = newConnection(jdbcConnection)
    try {
      val liquibase = new Liquibase(DefaultChangeLog, new ClassLoaderResourceAccessor(), liquibaseConnection)
      checkForChangeLogDir(liquibase)
      updateSchema(liquibase)
    } finally {
      closeConnection(liquibaseConnection)
    }
  }

  /**
    * Wraps a jdbc connection in the database with the appropriate liquibase connection.
    * As of 3.4.x, liquibase uses a custom connection for Hsql, Sybase, and Derby, although only Hsql is supported by
    * cromwell.
    *
    * @param jdbcConnection The liquibase connection.
    * @return
    */
  def newConnection(jdbcConnection: Connection): DatabaseConnection = {
    jdbcConnection.getMetaData.getDatabaseProductName match {
      case HsqlDatabaseProperties.PRODUCT_NAME => new HsqlConnection(jdbcConnection)
      case _ => new JdbcConnection(jdbcConnection)
    }
  }

  /**
    * Looks for evidence that the change log was previously applied using files in a directory under
    * src/main/migrations/. If found, throws an error stating that the database must be manually updated first to
    * remove the paths from the database. The sql returned in the error message updates the file paths, and also resets
    * the md5s to null such that they'll be recalculated on the next update.
    *
    * Good example:
    *
    *   mysql> select distinct filename from databasechangelog limit 3;
    *   +--------------------------------------+
    *   | filename                             |
    *   +--------------------------------------+
    *   | changesets/db_schema.xml             |
    *   | changesets/symbol_iteration_null.xml |
    *   | changesets/wdl_and_inputs.xml        |
    *   +--------------------------------------+
    *   3 rows in set (0.00 sec)
    *
    * Example that needs an update:
    *
    *   mysql> select distinct filename from databasechangelog limit 3;
    *   +---------------------------------------------------------------------------------------------------+
    *   | filename                                                                                          |
    *   +---------------------------------------------------------------------------------------------------+
    *   | /root/github.com/broadinstitute/cromwell/src/main/migrations/changesets/db_schema.xml             |
    *   | /root/github.com/broadinstitute/cromwell/src/main/migrations/changesets/symbol_iteration_null.xml |
    *   | /root/github.com/broadinstitute/cromwell/src/main/migrations/changesets/wdl_and_inputs.xml        |
    *   +---------------------------------------------------------------------------------------------------+
    *   3 rows in set (0.00 sec)
    *
    * @param liquibase The facade for interacting with liquibase.
    */
  def checkForChangeLogDir(liquibase: Liquibase): Unit = {
    val unexpectedChangeSets = liquibase.listUnexpectedChangeSets(DefaultContexts, DefaultLabelExpression).asScala
    val changeLogDirSets = unexpectedChangeSets.filter(_.getChangeLog.contains(ChangeLogDir))
    if (changeLogDirSets.nonEmpty) {
      throw new Error(
        s"""Backup and then run the following SQL command on your database before proceeding:
            |
            |    update DATABASECHANGELOG
            |    set MD5SUM = null,
            |      FILENAME = substr(FILENAME, instr(FILENAME, "$ChangeLogDir") + length("$ChangeLogDir"))
            |    where FILENAME like '%$ChangeLogDir%'
            |
            |Changesets were detected with path $ChangeLogDir:
            |    ${unexpectedChangeSets.mkString("\n    ")}
            |""".stripMargin)
    }
  }

  /**
    * Updates the liquibase database.
    *
    * @param liquibase The facade for interacting with liquibase.
    */
  def updateSchema(liquibase: Liquibase): Unit = {
    liquibase.update(DefaultContexts, DefaultLabelExpression)
  }

  /**
    * Converts a liquibase connection to a liquibase database.
    *
    * @param liquibaseConnection The liquibase connection.
    * @return The liquibase database.
    */
  def toDatabase(liquibaseConnection: DatabaseConnection): Database = {
    DatabaseFactory.getInstance().findCorrectDatabaseImplementation(liquibaseConnection)
  }

  /**
    * Compares a reference to a comparison liquibase database.
    *
    * @param referenceDatabase The reference liquibase database.
    * @param comparisonDatabase The comparison liquibase database.
    * @return The complete diff results.
    */
  def compare(referenceDatabase: Database, comparisonDatabase: Database): DiffResult = {
    DiffGeneratorFactory.getInstance().compare(referenceDatabase, comparisonDatabase, CompareControl.STANDARD)
  }

  /**
    * Compares a reference to a comparison JDBC connection.
    *
    * @param referenceJdbc The reference connection.
    * @param comparisonJdbc The comparison connection.
    * @return The complete diff results.
    */
  def compare(referenceJdbc: Connection, comparisonJdbc: Connection): DiffResult = {
    withConnection(referenceJdbc) { referenceLiquibase =>
      withConnection(comparisonJdbc) { comparisonLiquibase =>
        compare(toDatabase(referenceLiquibase), toDatabase(comparisonLiquibase))
      }
    }
  }

  /**
    * Provides a connection to a block of code, closing the connection afterwards.
    *
    * @param jdbcConnection The connection.
    * @param block The block to run.
    * @tparam T The return type of the block.
    * @return The result of running the block.
    */
  def withConnection[T](jdbcConnection: Connection)(block: DatabaseConnection => T): T = {
    val liquibaseConnection = newConnection(jdbcConnection)
    try {
      block(liquibaseConnection)
    } finally {
      closeConnection(liquibaseConnection)
    }
  }

  /**
    * Attempts to close a liquibase connection.
    *
    * @param connection The liquibase connection.
    */
  def closeConnection(connection: DatabaseConnection): Unit = {
    try {
      connection.close()
    } finally {
      /* ignore */
    }
  }
}

