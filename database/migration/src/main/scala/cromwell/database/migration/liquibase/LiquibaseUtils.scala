package cromwell.database.migration.liquibase

import java.sql.Connection

import liquibase.changelog.{ChangeLogParameters, ChangeSet, DatabaseChangeLog}
import liquibase.database.jvm.{HsqlConnection, JdbcConnection}
import liquibase.database.{Database, DatabaseConnection, DatabaseFactory, ObjectQuotingStrategy}
import liquibase.diff.compare.CompareControl
import liquibase.diff.{DiffGeneratorFactory, DiffResult}
import liquibase.parser.ChangeLogParserFactory
import liquibase.resource.ClassLoaderResourceAccessor
import liquibase.snapshot.{DatabaseSnapshot, SnapshotControl, SnapshotGeneratorFactory}
import liquibase.{Contexts, LabelExpression, Liquibase}
import org.hsqldb.persist.HsqlDatabaseProperties

import scala.collection.JavaConverters._

object LiquibaseUtils {
  // Paranoia: Create our own mutex. https://stackoverflow.com/questions/442564/avoid-synchronizedthis-in-java
  private val mutex = new Object
  private val DefaultContexts = new Contexts()
  private val DefaultLabelExpression = new LabelExpression()

  /**
    * Updates a liquibase schema to the latest version.
    *
    * @param settings The liquibase settings.
    * @param jdbcConnection A jdbc connection to the database.
    */
  def updateSchema(settings: LiquibaseSettings)(jdbcConnection: Connection): Unit = {
    mutex.synchronized {
      val liquibaseConnection = newConnection(jdbcConnection)
      try {
        val database = DatabaseFactory.getInstance.findCorrectDatabaseImplementation(liquibaseConnection)
        database.setDatabaseChangeLogLockTableName(settings.databaseChangeLogLockTableName.toUpperCase)
        database.setDatabaseChangeLogTableName(settings.databaseChangeLogTableName.toUpperCase)

        val liquibase = new Liquibase(settings.changeLogResourcePath, new ClassLoaderResourceAccessor(), database)
        updateSchema(liquibase)
      } finally {
        closeConnection(liquibaseConnection)
      }
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
  private def newConnection(jdbcConnection: Connection): DatabaseConnection = {
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
  private def updateSchema(liquibase: Liquibase): Unit = {
    liquibase.update(DefaultContexts, DefaultLabelExpression)
  }

  /**
    * Converts a liquibase connection to a liquibase database.
    *
    * @param liquibaseConnection The liquibase connection.
    * @return The liquibase database.
    */
  private def toDatabase(liquibaseConnection: DatabaseConnection): Database = {
    DatabaseFactory.getInstance().findCorrectDatabaseImplementation(liquibaseConnection)
  }

  /**
    * Compares a reference to a comparison liquibase database.
    *
    * @param referenceDatabase The reference liquibase database.
    * @param comparisonDatabase The comparison liquibase database.
    * @return The complete diff results.
    */
  private def compare(referenceDatabase: Database, comparisonDatabase: Database): DiffResult = {
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
    mutex.synchronized {
      withConnection(referenceJdbc) { referenceLiquibase =>
        withConnection(comparisonJdbc) { comparisonLiquibase =>
          val diffResult = compare(toDatabase(referenceLiquibase), toDatabase(comparisonLiquibase))
          block(diffResult)
        }
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
  private def withConnection[T](jdbcConnection: Connection)(block: DatabaseConnection => T): T = {
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
  private def closeConnection(connection: DatabaseConnection): Unit = {
    try {
      connection.close()
    } finally {
      /* ignore */
    }
  }

  /**
    * Returns the changelog for a liquibase setting.
    *
    * @param settings The liquibase settings.
    * @return The database changelog.
    */
  private def getChangeLog(settings: LiquibaseSettings): DatabaseChangeLog = {
    val changeLogFile: String = settings.changeLogResourcePath
    val resourceAccessor = new ClassLoaderResourceAccessor()
    val changeLogParameters = new ChangeLogParameters()
    val parser = ChangeLogParserFactory.getInstance.getParser(changeLogFile, resourceAccessor)
    val databaseChangeLog = parser.parse(changeLogFile, changeLogParameters, resourceAccessor)
    databaseChangeLog
  }

  /**
    * Returns the change sets for a liquibase setting.
    *
    * @param settings The liquibase settings.
    * @return The database change sets.
    */
  def getChangeSets(settings: LiquibaseSettings): Seq[ChangeSet] = {
    mutex.synchronized {
      getChangeLog(settings).getChangeSets.asScala
    }
  }

  /**
    * Returns a schema snapshot.
    *
    * @param jdbcConnection A jdbc connection to the database.
    * @return The database change sets.
    */
  def getSnapshot(jdbcConnection: Connection): DatabaseSnapshot = {
    mutex.synchronized {
      withConnection(jdbcConnection) { referenceLiquibase =>
        val database = toDatabase(referenceLiquibase)
        val objectQuotingStrategy = database.getObjectQuotingStrategy
        try {
          // Quote all objects for PostgreSQL
          database.setObjectQuotingStrategy(ObjectQuotingStrategy.QUOTE_ALL_OBJECTS)

          SnapshotGeneratorFactory.getInstance.createSnapshot(
            database.getDefaultSchema,
            database,
            new SnapshotControl(database)
          )
        } finally {
          database.setObjectQuotingStrategy(objectQuotingStrategy)
        }
      }
    }
  }
}
