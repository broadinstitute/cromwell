package cromwell.engine.db.slick

import java.sql.Connection

import _root_.slick.util.ConfigExtensionMethods._
import liquibase.Liquibase
import liquibase.database.DatabaseConnection
import liquibase.resource.{FileSystemResourceAccessor, ResourceAccessor}
import org.slf4j.LoggerFactory

import scala.concurrent.{ExecutionContext, Future}

class TestSlickDatabase(configPath: String) {

  private lazy val databaseConfig = DatabaseConfig.rootDatabaseConfig.getConfig(configPath)
  private lazy val log = LoggerFactory.getLogger(classOf[TestSlickDatabase])

  val dataAccessComponent: DataAccessComponent = new DataAccessComponent(databaseConfig.getString("slick.driver"))

  import dataAccessComponent.driver.api._

  /**
   * Check the database connection.
   * Can be run before operations that use slickDataAccess,
   * but creates a whole new connection pool to do so.
   */
  def isValidConnection: Future[Boolean] = {
    implicit val executionContext = ExecutionContext.global
    Future {
      log.debug("Opening test connection setup for " + configPath)
      Database.forConfig("", databaseConfig)
    } flatMap { database =>
      database.run(SimpleDBIO(_.connection.isValid(1))) recover {
        case ex =>
          log.error("Unable to connect to database under config: " + configPath, ex)
          false
      } andThen {
        case _ =>
          log.debug("Closing test connection setup for " + configPath)
          database.close()
      }
    }
  }

  lazy val slickDataAccess =
    if (this.databaseConfig == DatabaseConfig.databaseConfig)
      new SlickDataAccess() // Test the no-args constructor
    else
      new SlickDataAccess(databaseConfig, dataAccessComponent)

  def useLiquibase = databaseConfig.hasPath("liquibase")

  def setupLiquibase(): Unit = {
    val liquibaseContexts = "test"
    val liquibaseChangeLog = databaseConfig.getString("liquibase.changelog")
    val liquibaseConnectionClass = databaseConfig.getStringOr(
      "liquibase.connection", "liquibase.database.jvm.JdbcConnection")
    val liquibaseDropAll = databaseConfig.getBooleanOr("liquibase.dropall")

    val connectionClass = Class.forName(liquibaseConnectionClass)
    val connectionConstructor = connectionClass.getConstructor(classOf[Connection])
    val jdbcConnection = slickDataAccess.database.source.createConnection()
    val liquibaseConnection: DatabaseConnection =
      connectionConstructor.newInstance(jdbcConnection).asInstanceOf[DatabaseConnection]
    try {
      val resourceAccessor: ResourceAccessor = new FileSystemResourceAccessor()
      val liquibase = new Liquibase(liquibaseChangeLog, resourceAccessor, liquibaseConnection)
      if (liquibaseDropAll)
        liquibase.dropAll()
      liquibase.update(liquibaseContexts)
    } finally {
      liquibaseConnection.close()
      jdbcConnection.close()
    }
  }

}
