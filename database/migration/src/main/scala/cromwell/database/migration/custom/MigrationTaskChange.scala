package cromwell.database.migration.custom

import com.typesafe.config.ConfigFactory
import com.typesafe.scalalogging.LazyLogging
import liquibase.change.custom.CustomTaskChange
import liquibase.database.Database
import liquibase.database.jvm.JdbcConnection
import liquibase.exception.{CustomChangeException, ValidationErrors}
import liquibase.resource.ResourceAccessor

/**
  * Provides a default implementation of a liquibase custom task change.
  */
trait MigrationTaskChange extends CustomTaskChange with LazyLogging {
  lazy val config = ConfigFactory.load

  /** @return name of the migration, defaulting to the class name */
  def migrationName: String = getClass.getSimpleName

  /**
    * Performs the migration.
    *
    * @param connection the connection to the database
    */
  def migrate(connection: JdbcConnection): Unit

  override def execute(database: Database): Unit = {
    try {
      val dbConn = database.getConnection.asInstanceOf[JdbcConnection]
      val autoCommit = dbConn.getAutoCommit
      dbConn.setAutoCommit(false)
      migrate(dbConn)
      dbConn.setAutoCommit(autoCommit)
    } catch {
      case customChangeException: CustomChangeException => throw customChangeException
      case exception: Exception =>
        throw new CustomChangeException(s"Could not apply migration script for $migrationName", exception)
    }
  }

  override def setUp() = {}

  override def getConfirmationMessage = s"$migrationName complete."

  override def validate(database: Database) = new ValidationErrors

  override def setFileOpener(resourceAccessor: ResourceAccessor) = {}
}
