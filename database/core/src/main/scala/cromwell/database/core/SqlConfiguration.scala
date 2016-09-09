package cromwell.database.core

import com.typesafe.config.ConfigFactory
import lenthall.config.ScalaConfig._

object SqlConfiguration {
  lazy val rootConfig = ConfigFactory.load()
  private lazy val rootDatabaseConfig = rootConfig.getConfig("database")
  private lazy val databaseConfigName = rootDatabaseConfig.getStringOption("config")
  lazy val defaultDatabaseConfig = databaseConfigName.map(getDatabaseConfig).getOrElse(rootDatabaseConfig)

  def getDatabaseConfig(path: String) = rootDatabaseConfig.getConfig(path)
}
