package cromwell.services.database

import java.sql.Connection

import better.files._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.database.migration.liquibase.LiquibaseUtils
import cromwell.database.slick.SlickDatabase
import cromwell.services.{EngineServicesStore, MetadataServicesStore}
import liquibase.snapshot.DatabaseSnapshot
import liquibase.structure.core.Index
import slick.jdbc.JdbcProfile
import slick.jdbc.meta.{MIndexInfo, MPrimaryKey}
import cromwell.services.ServicesStore.EnhancedSqlDatabase

import scala.concurrent.Await
import scala.concurrent.duration.Duration

object DatabaseTestKit {

  /**
    * Lends a connection to a block of code.
    *
    * @param profile  The slick jdbc profile for accessing the database.
    * @param database The database to use for the connection.
    * @param block    The block of code to run over the connection.
    * @tparam Profile The slick jdbc profile for accessing the database.
    * @tparam A       The return type of the block.
    * @return The return value of the block.
    */
  def withConnection[Profile <: JdbcProfile, A](profile: Profile, database: Profile#Backend#Database)
                                               (block: Connection => A): A = {
    /*
     TODO: Should this withConnection() method have a (implicit?) timeout parameter, that it passes on to Await.result?
     If we run completely asynchronously, nest calls to withConnection, and then call flatMap, the outer connection may
     already be closed before an inner block finishes running.
     */
    Await.result(database.run(profile.api.SimpleDBIO(context => block(context.connection))), Duration.Inf)
  }

  /**
    * Lends two connections to a block of code.
    *
    * @param profile1  The slick jdbc profile for accessing the first database.
    * @param database1 The database to use for the first connection.
    * @param profile2  The slick jdbc profile for accessing the second database.
    * @param database2 The database to use for the second connection.
    * @param block     The block of code to run over the first and second connections.
    * @tparam Profile1 The slick jdbc profile for accessing the first database.
    * @tparam Profile2 The slick jdbc profile for accessing the second database.
    * @tparam A        The return type of the block.
    * @return The return value of the block.
    */
  def withConnections[Profile1 <: JdbcProfile, Profile2 <: JdbcProfile, A]
  (profile1: Profile1, database1: Profile1#Backend#Database, profile2: Profile2, database2: Profile2#Backend#Database)
  (block: (Connection, Connection) => A): A = {
    withConnection(profile1, database1) { connection1 =>
      withConnection(profile2, database2) { connection2 =>
        block(connection1, connection2)
      }
    }
  }

  /**
    * Creates a new in memory HSQLDB that should be closed after use.
    */
  def inMemoryDatabase[A <: SlickDatabase](databaseType: CromwellDatabaseType[A], schemaManager: SchemaManager): A = {
    val databaseConfig = ConfigFactory.parseString(
      s"""|db.url = "jdbc:hsqldb:mem:$${uniqueSchema};shutdown=false;hsqldb.tx=mvcc"
          |db.driver = "org.hsqldb.jdbcDriver"
          |db.connectionTimeout = 3000
          |profile = "slick.jdbc.HsqldbProfile$$"
          |liquibase.updateSchema = false
          |""".stripMargin)
    val database = databaseType.newDatabase(databaseConfig)
    schemaManager match {
      case SlickSchemaManager => SlickDatabase.createSchema(database)
      case LiquibaseSchemaManager => runLiquibase(database, databaseType)
    }
    database
  }

  /**
    * Opens an initialized database.
    */
  def initializedDatabaseFromConfig[A <: SlickDatabase](databaseType: CromwellDatabaseType[A],
                                                        databaseConfig: Config): A = {
    val database = databaseType.newDatabase(databaseConfig)
    database.initialized(databaseType.liquibaseSettings)
  }

  /**
    * Opens an initialized database.
    */
  def initializedDatabaseFromSystem[A <: SlickDatabase](databaseType: CromwellDatabaseType[A],
                                                        databaseSystem: DatabaseSystem): A = {
    val databaseConfig = ConfigFactory.load.getConfig(databaseSystem.configPath)
    initializedDatabaseFromConfig(databaseType, databaseConfig)
  }

  /**
    * Run liquibase on a open database.
    */
  def runLiquibase(database: SlickDatabase, databaseType: CromwellDatabaseType[_]): Unit = {
    val settings = databaseType match {
      case EngineDatabaseType => EngineServicesStore.EngineLiquibaseSettings
      case MetadataDatabaseType => MetadataServicesStore.MetadataLiquibaseSettings
    }
    database withConnection LiquibaseUtils.updateSchema(settings)
  }

  /**
    * Returns a Liquibase snapshot of an open Slick database.
    */
  def liquibaseSnapshot(database: SlickDatabase): DatabaseSnapshot = {
    withConnection(database.dataAccess.driver, database.database)(LiquibaseUtils.getSnapshot)
  }

  /**
    * Returns a Liquibase snapshot of an in memory HSQLDB.
    */
  def inMemorySnapshot[A <: SlickDatabase](databaseType: CromwellDatabaseType[A],
                       schemaManager: SchemaManager): DatabaseSnapshot = {
    var snapshot: DatabaseSnapshot = null
    for {
      database <- inMemoryDatabase(databaseType, schemaManager).autoClosed
    } {
      snapshot = liquibaseSnapshot(database)
    }
    snapshot
  }

  /**
    * Returns true if the primary key was auto generated by the database.
    */
  def isGenerated(primaryKey: MPrimaryKey): Boolean = {
    isGenerated(primaryKey.pkName.get)
  }

  /**
    * Returns true if the index was auto generated by the database.
    */
  def isGenerated(index: MIndexInfo): Boolean = {
    isGenerated(index.indexName.get)
  }

  /**
    * Returns true if the index was auto generated by the database.
    */
  def isGenerated(index: Index): Boolean = {
    isGenerated(index.getName)
  }

  private def isGenerated(name: String): Boolean = {
    name.startsWith("SYS_")
  }
}
