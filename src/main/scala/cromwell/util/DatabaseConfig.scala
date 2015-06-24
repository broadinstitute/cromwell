package cromwell.util

import com.typesafe.config.ConfigFactory
import slick.util.ConfigExtensionMethods._

object DatabaseConfig {
  private val config = ConfigFactory.load()
  private val rootDatabaseConfig = config.getConfig("database")
  private val databaseConfigName = rootDatabaseConfig.getStringOr("config", null)
  lazy val databaseConfig = if (databaseConfigName == null) rootDatabaseConfig else rootDatabaseConfig.getConfig(databaseConfigName)
  lazy val slickDriver = databaseConfig.getString("slick.driver")
  lazy val liquibaseSetup = databaseConfig.hasPath("liquibase")
  lazy val liquibaseChangeLog = databaseConfig.getString("liquibase.changelog")
  lazy val liquibaseConnection = databaseConfig.getStringOr("liquibase.connection", "liquibase.database.jvm.JdbcConnection")
  lazy val liquibaseDropAll = databaseConfig.getBooleanOr("liquibase.dropall")
}
