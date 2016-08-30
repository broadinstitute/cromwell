package cromwell.database.migration.restart.table

import liquibase.change.custom.CustomTaskChange
import liquibase.database.Database
import liquibase.database.jvm.JdbcConnection
import liquibase.exception.ValidationErrors
import liquibase.resource.ResourceAccessor

abstract class AbstractRestartMigration extends CustomTaskChange {

  protected def description: String

  protected def doMigration(connection: JdbcConnection): Unit

  override def execute(database: Database): Unit = {
    val connection = database.getConnection.asInstanceOf[JdbcConnection]
    connection.setAutoCommit(false)
    doMigration(connection)
    connection.commit()
  }

  override def setUp(): Unit = ()

  override def getConfirmationMessage: String = s"$description migration complete."

  override def validate(database: Database): ValidationErrors = new ValidationErrors

  override def setFileOpener(resourceAccessor: ResourceAccessor): Unit = ()
}
