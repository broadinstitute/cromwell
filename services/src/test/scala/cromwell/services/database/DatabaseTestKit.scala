package cromwell.services.database

import java.net.URLEncoder
import java.sql.Connection

import better.files._
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
    val databaseConfig = getConfig(databaseSystem)
    initializedDatabaseFromConfig(databaseType, databaseConfig)
  }

  /**
    * Opens a database connection without any liquibase being performed.
    */
  def schemalessDatabaseFromSystem(databaseSystem: DatabaseSystem): SchemalessSlickDatabase with TestSlickDatabase = {
    val databaseConfig = getConfig(databaseSystem)
    new SchemalessSlickDatabase(databaseConfig) with TestSlickDatabase
  }

  private var configCache: Map[NetworkDatabaseSystem, Config] = Map.empty
  private val configCacheMutex = new Object

  /**
    * Returns a config for a DatabaseSystem.
    */
  private def getConfig(databaseSystem: DatabaseSystem): Config = {
    databaseSystem match {
      case HsqldbDatabaseSystem => hsqldbDatabaseConfig
      case networkDatabaseSystem: NetworkDatabaseSystem =>
        configCacheMutex synchronized {
          configCache.get(networkDatabaseSystem) match {
            case Some(config) => config
            case None =>
              val config = getConfig(networkDatabaseSystem)
              configCache += networkDatabaseSystem -> config
              config
          }
        }
    }
  }

  private case class DatabaseSystemSettings(environmentKey: String, defaultPort: Int, dockerTag: String)

  /**
    * Returns the network settings for a database system.
    */
  private def getDatabaseSystemSettings(networkDatabaseSystem: NetworkDatabaseSystem): DatabaseSystemSettings = {
    networkDatabaseSystem match {
      // The below list of docker tags should be synced with the tags under "BUILD_TYPE=dbms" in .travis.yml
      case MariadbEarliestDatabaseSystem => DatabaseSystemSettings("MARIADB", 23306, "5.5")
      case MariadbLatestDatabaseSystem => DatabaseSystemSettings("MARIADB_LATEST", 33306, "latest")
      case MysqlEarliestDatabaseSystem => DatabaseSystemSettings("MYSQL", 3306, "5.6")
      case MysqlLatestDatabaseSystem => DatabaseSystemSettings("MYSQL_LATEST", 13306, "latest")
      case PostgresqlEarliestDatabaseSystem => DatabaseSystemSettings("POSTGRESQL", 5432, "9.6")
      case PostgresqlLatestDatabaseSystem => DatabaseSystemSettings("POSTGRESQL_LATEST", 15432, "latest")
      // The list above of docker tags should be synced with the tags under "BUILD_TYPE=dbms" in .travis.yml
    }
  }

  /**
    * Returns a config for a NetworkDatabaseSystem.
    */
  private def getConfig(networkDatabaseSystem: NetworkDatabaseSystem): Config = {
    val databaseSystemSettings = getDatabaseSystemSettings(networkDatabaseSystem)
    val systemName = networkDatabaseSystem.name
    val environmentKey = databaseSystemSettings.environmentKey
    val jdbcPortDefault = databaseSystemSettings.defaultPort
    val dockerTag = databaseSystemSettings.dockerTag

    val jdbcUsername = "cromwell"
    val jdbcPassword = "test"
    val jdbcSchema = "cromwell_test"
    val jdbcHostname: String = sys.env.getOrElse(s"CROMWELL_BUILD_${environmentKey}_HOSTNAME", "localhost")
    val jdbcPort = sys.env.get(s"CROMWELL_BUILD_${environmentKey}_PORT").map(_.toInt).getOrElse(jdbcPortDefault)

    def makeJdbcUrl(dbms: String, queryParams: Map[String, String]): String = {
      s"jdbc:$dbms://$jdbcHostname:$jdbcPort/$jdbcSchema?" +
        queryParams.map({ case (name, value) => queryEncode(name) + "=" + queryEncode(value) }).mkString("&")
    }

    val (dockerHelp, resetHelp, slickProfile, jdbcDriver, jdbcUrl) = networkDatabaseSystem.platform match {
      case HsqldbDatabasePlatform => throw new UnsupportedOperationException
      case MariadbDatabasePlatform => (
        s"""|docker run \\
            |  --detach --name cromwell_database_$jdbcPort \\
            |  --env MYSQL_ROOT_PASSWORD=private \\
            |  --env MYSQL_USER=$jdbcUsername \\
            |  --env MYSQL_PASSWORD=$jdbcPassword \\
            |  --env MYSQL_DATABASE=$jdbcSchema \\
            |  --publish $jdbcPort:3306 \\
            |  --volume $${PWD}/src/ci/docker-compose/mariadb-conf.d:/etc/mysql/conf.d \\
            |  mariadb:$dockerTag
            |""".stripMargin.trim,
        s"""|mysql \\
            |  --protocol=tcp --host=$jdbcHostname --port=$jdbcPort \\
            |  --user=$jdbcUsername --password=$jdbcPassword \\
            |  --execute='DROP DATABASE IF EXISTS $jdbcSchema; CREATE DATABASE $jdbcSchema;'
            |""".stripMargin.trim,
        "slick.jdbc.MySQLProfile$",
        "org.mariadb.jdbc.Driver",
        makeJdbcUrl("mysql", Map("rewriteBatchedStatements" -> "true")),
      )
      case MysqlDatabasePlatform => (
        s"""|docker run \\
            |  --detach --name cromwell_database_$jdbcPort \\
            |  --env MYSQL_ROOT_PASSWORD=private \\
            |  --env MYSQL_USER=$jdbcUsername \\
            |  --env MYSQL_PASSWORD=$jdbcPassword \\
            |  --env MYSQL_DATABASE=$jdbcSchema \\
            |  --publish $jdbcPort:3306 \\
            |  --volume $${PWD}/src/ci/docker-compose/mysql-conf.d:/etc/mysql/conf.d \\
            |  mysql:$dockerTag
            |""".stripMargin.trim,
        s"""|mysql \\
            |  --protocol=tcp --host=$jdbcHostname --port=$jdbcPort \\
            |  --user=$jdbcUsername --password=$jdbcPassword \\
            |  --execute='DROP DATABASE IF EXISTS $jdbcSchema; CREATE DATABASE $jdbcSchema;'
            |""".stripMargin.trim,
        "slick.jdbc.MySQLProfile$",
        "com.mysql.cj.jdbc.Driver",
        makeJdbcUrl("mysql", Map(
          "rewriteBatchedStatements" -> "true",
          "useSSL" -> "false",
          "allowPublicKeyRetrieval" -> "true",
          "serverTimezone" -> "UTC",
          "useInformationSchema" -> "true",
        )),
      )
      case PostgresqlDatabasePlatform => (
        s"""|docker run \\
            |  --detach --name cromwell_database_$jdbcPort \\
            |  --env POSTGRES_USER=$jdbcUsername \\
            |  --env POSTGRES_PASSWORD=$jdbcPassword \\
            |  --env POSTGRES_DB=$jdbcSchema \\
            |  --publish $jdbcPort:5432 \\
            |  --volume $${PWD}/src/ci/docker-compose/postgresql-initdb.d:/docker-entrypoint-initdb.d \\
            |  postgres:$dockerTag
            |""".stripMargin.trim,
        s"""|PGPASSWORD=$jdbcPassword psql \\
            |  --host=localhost --port=5432 --username=$jdbcUsername \\
            |  postgres <<< 'drop database if exists $jdbcSchema; create database $jdbcSchema;'
            |""".stripMargin.trim,
        "slick.jdbc.PostgresProfile$",
        "org.postgresql.Driver",
        makeJdbcUrl("postgresql", Map("reWriteBatchedInserts" -> "true")),
      )
    }

    logger.info(
      s"""|Run an example $systemName via docker using:
          |$dockerHelp""".stripMargin)
    logger.info(
      s"""|The schema will be initialized when the docker container starts. If needed reset the schema using:
          |$resetHelp""".stripMargin)

    ConfigFactory.parseString(
      s"""|profile = "$slickProfile"
          |db {
          |  driver = "$jdbcDriver"
          |  url = "$jdbcUrl"
          |  user = "$jdbcUsername"
          |  password = "$jdbcPassword"
          |  connectionTimeout = 5000
          |}
          |""".stripMargin
    )
  }

  /**
    * Encode strings for URLs.
    */
  private def queryEncode(string: String): String = URLEncoder.encode(string, "UTF-8")

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
