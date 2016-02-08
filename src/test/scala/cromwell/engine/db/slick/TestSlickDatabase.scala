package cromwell.engine.db.slick

import org.slf4j.LoggerFactory
import slick.backend.DatabaseConfig
import slick.driver.JdbcProfile

import scala.concurrent.{ExecutionContext, Future}

class TestSlickDatabase(configPath: String) {

  private lazy val databaseConfig = SlickDataAccess.getDatabaseConfig(configPath)
  private lazy val log = LoggerFactory.getLogger(classOf[TestSlickDatabase])

  private lazy val slickConfig = DatabaseConfig.forConfig[JdbcProfile]("", databaseConfig)

  // NOTE: Using the import below for isValidConnection, but maybe not the lazy instance if the check fails.
  lazy val dataAccessComponent: DataAccessComponent = new DataAccessComponent(slickConfig.driver)

  import dataAccessComponent.driver.api._

  /**
   * Check the database connection.
   *
   * This check only produces warnings, not errors. The primary use case is checking if _optional_ tests can be run
   * against a database configuration.
   *
   * Can be run before operations that use slickDataAccess, but creates a whole new connection pool to do so.
   */
  def isValidConnection: Future[Boolean] = {
    implicit val executionContext = ExecutionContext.global
    Future {
      log.debug("Opening test connection setup for " + configPath)
      slickConfig.db
    } flatMap { database =>
      database.run(SimpleDBIO(_.connection.isValid(1))) andThen {
        case _ =>
          log.debug("Closing test connection setup for " + configPath)
          database.close()
      }
    } recover {
      case ex =>
        log.warn("Unable to connect to database under config: " + configPath, ex)
        false
    }
  }

  lazy val slickDataAccess =
    if (this.databaseConfig == SlickDataAccess.defaultDatabaseConfig)
      new SlickDataAccess() // Test the no-args constructor
    else
      new SlickDataAccess(databaseConfig)
}
