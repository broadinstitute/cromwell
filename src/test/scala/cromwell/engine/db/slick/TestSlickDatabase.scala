package cromwell.engine.db.slick

import java.sql.Connection

import _root_.slick.util.ConfigExtensionMethods._
import com.typesafe.config.Config
import liquibase.Liquibase
import liquibase.database.DatabaseConnection
import liquibase.resource.{FileSystemResourceAccessor, ResourceAccessor}

class TestSlickDatabase(databaseConfig: Config) {

  lazy val slickDataAccess = new SlickDataAccess(databaseConfig)

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
