package cromwell.services.database

import java.sql.Connection

import better.files._
import com.dimafeng.testcontainers.{Container, JdbcDatabaseContainer, MariaDBContainer, MySQLContainer, PostgreSQLContainer}
import com.typesafe.config.{Config, ConfigFactory}
import com.typesafe.scalalogging.StrictLogging
import cromwell.database.migration.liquibase.LiquibaseUtils
import cromwell.database.slick.SlickDatabase
import cromwell.services.ServicesStore.EnhancedSqlDatabase
import cromwell.services.{EngineServicesStore, MetadataServicesStore}
import liquibase.snapshot.DatabaseSnapshot
import liquibase.structure.core.Index
import slick.jdbc.JdbcProfile
import slick.jdbc.meta.{MIndexInfo, MPrimaryKey}

import scala.concurrent.Await
import scala.concurrent.duration.Duration

object DatabaseTestKit extends StrictLogging {

  private lazy val hsqldbDatabaseConfig = ConfigFactory.load().getConfig("database")

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
    val databaseConfig = ConfigFactory.parseString("liquibase.updateSchema = false") withFallback hsqldbDatabaseConfig
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
                                                        databaseConfig: Config): A with TestSlickDatabase = {
    val database = databaseType.newDatabase(databaseConfig)
    database.initialized(databaseType.liquibaseSettings)
  }

  /**
    * Opens an initialized database.
    */
  def initializedDatabaseFromSystem[A <: SlickDatabase](databaseType: CromwellDatabaseType[A],
                                                        databaseSystem: DatabaseSystem): A with TestSlickDatabase = {
    initializeDatabaseByContainerOptTypeAndSystem(None, databaseType, databaseSystem)
  }

  /**
   * Opens a database connection without any liquibase being performed.
   */
  def schemalessDatabaseFromContainerOptAndSystem(containerOpt: Option[Container], databaseSystem: DatabaseSystem): SchemalessSlickDatabase with TestSlickDatabase = {
    containerOpt match {
      case None => new SchemalessSlickDatabase(getConfig(databaseSystem, dbContainerOpt = None)) with TestSlickDatabase
      case Some(cont) if cont.isInstanceOf[JdbcDatabaseContainer] =>
        new SchemalessSlickDatabase(getConfig(databaseSystem, Option(cont.asInstanceOf[JdbcDatabaseContainer]))) with TestSlickDatabase
      case Some(_) => throw new RuntimeException("ERROR: container is not a JdbcDatabaseContainer.")
    }
  }

  def getDatabaseTestContainer(databaseSystem: DatabaseSystem): Option[Container] = {
    databaseSystem match {
      case HsqldbDatabaseSystem => None
      case networkDbSystem: NetworkDatabaseSystem =>
        networkDbSystem.platform match {
          case MariadbDatabasePlatform =>
            Option(MariaDBContainer(
              dockerImageName = s"mariadb:${networkDbSystem.dockerImageVersion}",
              dbName = "cromwell_test",
              dbUsername = "cromwell",
              dbPassword = "test"
            ))
          case MysqlDatabasePlatform =>
            Option(MySQLContainer(
              mysqlImageVersion = s"mysql:${networkDbSystem.dockerImageVersion}",
              databaseName = "cromwell_test",
              username = "cromwell",
              password = "test"))
          case PostgresqlDatabasePlatform =>
            Option(PostgreSQLContainer(
              dockerImageNameOverride =  s"postgres:${networkDbSystem.dockerImageVersion}",
              databaseName = "cromwell_test",
              username = "cromwell",
              password = "test"))
          case _ => None
        }
    }
  }

  def initializeDatabaseByContainerOptTypeAndSystem[A <: SlickDatabase](containerOpt: Option[Container],
                                                                        databaseType: CromwellDatabaseType[A],
                                                                        databaseSystem: DatabaseSystem): A with TestSlickDatabase = {
    containerOpt match {
      case None => initializedDatabaseFromConfig(databaseType, getConfig(databaseSystem, None))
      case Some(cont) if cont.isInstanceOf[JdbcDatabaseContainer] =>
        initializedDatabaseFromConfig(databaseType, getConfig(databaseSystem, Option(cont.asInstanceOf[JdbcDatabaseContainer])))
      case Some(_) => throw new RuntimeException("ERROR: container is not a JdbcDatabaseContainer.")
    }
  }

  def getConfig(databaseSystem: DatabaseSystem, dbContainerOpt: Option[JdbcDatabaseContainer] = None): Config = {
    dbContainerOpt match {
      case None if databaseSystem == HsqldbDatabaseSystem => hsqldbDatabaseConfig
      case None => throw new RuntimeException("ERROR: dbContainer option must be passed into `getConfig` method for all databases except HSQLDB.")
      case Some(dbContainer) =>
        val (slickProfile, jdbcDriver) = databaseSystem.platform match {
          case HsqldbDatabasePlatform => throw new RuntimeException("ERROR: dbContainer option cannot be defined for HSQLDB database.")
          case MariadbDatabasePlatform => ("slick.jdbc.MySQLProfile$", "org.mariadb.jdbc.Driver")
          case MysqlDatabasePlatform => ("slick.jdbc.MySQLProfile$", "com.mysql.cj.jdbc.Driver")
          case PostgresqlDatabasePlatform => ("slick.jdbc.PostgresProfile$", "org.postgresql.Driver")
        }
        ConfigFactory.parseString(
          s"""|profile = "$slickProfile"
              |db {
              |  driver = "$jdbcDriver"
              |  url = "${dbContainer.jdbcUrl}"
              |  user = "${dbContainer.username}"
              |  password = "${dbContainer.password}"
              |  connectionTimeout = 5000
              |}
              |""".stripMargin
        )
    }
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
    * Returns the database connection metadata for an open Slick database.
    */
  def connectionMetadata(database: SlickDatabase): ConnectionMetadata = {
    withConnection(database.dataAccess.driver, database.database) {
      connection =>
        val metadata = connection.getMetaData
        ConnectionMetadata(
          databaseProductName = metadata.getDatabaseProductName,
          databaseProductVersion = metadata.getDatabaseProductVersion,
          databaseMajorVersion = metadata.getDatabaseMajorVersion,
          databaseMinorVersion = metadata.getDatabaseMinorVersion,
          driverName = metadata.getDriverName,
          driverVersion = metadata.getDriverVersion,
          driverMajorVersion = metadata.getDriverMajorVersion,
          driverMinorVersion = metadata.getDriverMinorVersion,
          jdbcMajorVersion = metadata.getJDBCMajorVersion,
          jdbcMinorVersion = metadata.getJDBCMinorVersion,
        )
    }
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
