package cromwell.engine.db.slick

import com.typesafe.config.ConfigFactory
import slick.util.ConfigExtensionMethods._

object DatabaseConfig {
  private val config = ConfigFactory.load()
  lazy val rootDatabaseConfig = config.getConfig("database")
  private val databaseConfigName = rootDatabaseConfig.getStringOpt("config")
  lazy val databaseConfig = databaseConfigName.map(rootDatabaseConfig.getConfig).getOrElse(rootDatabaseConfig)
}
