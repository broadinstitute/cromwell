package cromwell.database.sql

import java.sql.Connection
import java.util.UUID
import lenthall.config.ScalaConfig._
import com.typesafe.config.{Config, ConfigFactory, ConfigValueFactory}

trait SqlDatabase extends AutoCloseable
  with OldeWorldeSqlDatabase
  with MetadataSqlDatabase
  with WorkflowStoreSqlDatabase
  with BackendKVStoreSqlDatabase
  with JobStoreSqlDatabase
  with CallCachingStore {

  val urlKey: String
  val originalDatabaseConfig: Config
  lazy val databaseConfig = SqlDatabase.withUniqueSchema(originalDatabaseConfig, urlKey)

  def withConnection[A](block: Connection => A): A
}

object SqlDatabase {
  lazy val rootConfig = ConfigFactory.load()
  private lazy val rootDatabaseConfig = rootConfig.getConfig("database")
  private lazy val databaseConfigName = rootDatabaseConfig.getStringOption("config")
  lazy val defaultDatabaseConfig = databaseConfigName.map(getDatabaseConfig).getOrElse(rootDatabaseConfig)

  def getDatabaseConfig(path: String) = rootDatabaseConfig.getConfig(path)

  /**
    * Modifies config.getString("url") to return a unique schema, if the original url contains the text
    * "\${uniqueSchema}".
    *
    * This allows each instance of a database object to use a clean, and different, in memory database.
    *
    * @return Config with \${uniqueSchema} in url replaced with a unique string.
    */
  def withUniqueSchema(config: Config, urlKey: String): Config = {
    val urlValue = config.getString(urlKey)
    if (urlValue.contains("${uniqueSchema}")) {
      // Config wasn't updating with a simple withValue/withFallback.
      // So instead, do a bit of extra work to insert the generated schema name in the url.
      val schema = UUID.randomUUID().toString
      val newUrl = urlValue.replaceAll("""\$\{uniqueSchema\}""", schema)
      val origin = urlKey + " with uniqueSchema=" + schema
      val urlConfigValue = ConfigValueFactory.fromAnyRef(newUrl, origin)
      val urlConfig = ConfigFactory.empty(origin).withValue(urlKey, urlConfigValue)
      urlConfig.withFallback(config)
    } else {
      config
    }
  }
}
