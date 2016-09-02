package cromwell.database.migration.liquibase

import java.sql.Connection

import liquibase.database.jvm.{HsqlConnection, JdbcConnection}
import liquibase.database.{Database, DatabaseConnection, DatabaseFactory}
import liquibase.diff.compare.CompareControl
import liquibase.diff.{DiffGeneratorFactory, DiffResult}
import liquibase.resource.ClassLoaderResourceAccessor
import liquibase.{Contexts, LabelExpression, Liquibase}
import org.hsqldb.persist.HsqlDatabaseProperties

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
    * @param referenceJdbc  The reference connection.
    * @param comparisonJdbc The comparison connection.
    * @param block          Block of code to run before closing the connections.
    * @return The complete diff results.
    */
  def compare[T](referenceJdbc: Connection, comparisonJdbc: Connection)(block: DiffResult => T): T = {
    withConnection(referenceJdbc) { referenceLiquibase =>
      withConnection(comparisonJdbc) { comparisonLiquibase =>
        val diffResult = compare(toDatabase(referenceLiquibase), toDatabase(comparisonLiquibase))
        block(diffResult)
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
