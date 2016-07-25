package cromwell.database.slick

import java.util.UUID
import java.util.concurrent.{ExecutorService, Executors}

import com.typesafe.config.{Config, ConfigFactory, ConfigValueFactory}
import cromwell.database.slick.tables.DataAccessComponent
import cromwell.database.sql.SqlDatabase
import lenthall.config.ScalaConfig._
import org.slf4j.LoggerFactory
import slick.backend.DatabaseConfig
import slick.driver.JdbcProfile

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, Future}

object SlickDatabase {
  lazy val rootConfig = ConfigFactory.load()
  private lazy val rootDatabaseConfig = rootConfig.getConfig("database")
  private lazy val databaseConfigName = rootDatabaseConfig.getStringOption("config")
  lazy val defaultDatabaseConfig = databaseConfigName.map(getDatabaseConfig).getOrElse(rootDatabaseConfig)

  def getDatabaseConfig(path: String) = rootDatabaseConfig.getConfig(path)

  implicit class ConfigWithUniqueSchema(val config: Config) extends AnyVal {
    /**
      * Returns either the "url" or "properties.url"
      */
    def urlKey = if (config.hasPath("db.url")) "db.url" else "db.properties.url"

    /**
      * Returns the value of either the "url" or "properties.url"
      */
    def urlValue = config.getString(urlKey)

    /**
      * Modifies config.getString("url") to return a unique schema, if the original url contains the text
      * "\${slick.uniqueSchema}".
      *
      * This allows each instance of a SlickDataAccess object to use a clean, and different, in memory database.
      *
      * @return Config with \${slick.uniqueSchema} in url replaced with a unique string.
      */
    def withUniqueSchema: Config = {
      if (urlValue.contains("${slick.uniqueSchema}")) {
        // Config wasn't updating with a simple withValue/withFallback.
        // So instead, do a bit of extra work to insert the generated schema name in the url.
        val schema = UUID.randomUUID().toString
        val newUrl = urlValue.replaceAll("""\$\{slick\.uniqueSchema\}""", schema)
        val origin = urlKey + " with slick.uniqueSchema=" + schema
        val urlConfigValue = ConfigValueFactory.fromAnyRef(newUrl, origin)
        val urlConfig = ConfigFactory.empty(origin).withValue(urlKey, urlConfigValue)
        urlConfig.withFallback(config)
      } else {
        config
      }
    }
  }

  lazy val log = LoggerFactory.getLogger("cromwell.db.slick")
}

/**
  * Data Access implementation using Slick.
  *
  * NOTE: the uses of .head below will cause an exception to be thrown
  * if the list is empty.  In every use case as of the writing of this comment,
  * those exceptions would have been wrapped in a failed Future and returned.
  */
class SlickDatabase(databaseConfig: Config) extends SqlDatabase
  with OldeWorldeSlickDatabase
  with MetadataSlickDatabase
  with WorkflowStoreSlickDatabase
  with JobStoreSlickDatabase{

  import SlickDatabase._

  def this() = this(SlickDatabase.defaultDatabaseConfig)

  private val configWithUniqueSchema = this.databaseConfig.withUniqueSchema

  val slickConfig = DatabaseConfig.forConfig[JdbcProfile]("", configWithUniqueSchema)
  val dataAccess = new DataAccessComponent(slickConfig.driver)

  // Allows creation of a Database, plus implicits for running transactions
  import dataAccess.driver.api._

  // NOTE: if you want to refactor database is inner-class type: this.dataAccess.driver.backend.DatabaseFactory
  val database = slickConfig.db

  // Possibly create the database
  {
    import SlickDatabase._
    log.info(s"Running with database ${configWithUniqueSchema.urlKey} = ${configWithUniqueSchema.urlValue}")
    // NOTE: Slick 3.0.0 schema creation, Clobs, and MySQL don't mix:  https://github.com/slick/slick/issues/637
    //
    // Not really an issue, since externally run liquibase is standard way of installing / upgrading MySQL.
    //
    // Also, creating the unique key on UUID stored as a VARCHAR requires setting the length to O.Length(36) or (100)
    // for MySQL schema gen to avoid:
    //   com.mysql.jdbc.exceptions.jdbc4.MySQLSyntaxErrorException: BLOB/TEXT column 'WORKFLOW_EXECUTION_UUID'
    //   used in key specification without a key length
    //
    // Perhaps we'll use a more optimized data type for UUID's bytes in the future, as a FK, instead auto-inc cols
    //
    // The value `${slick.uniqueSchema}` may be used in the url, in combination with `slick.createSchema = true`, to
    // generate unique schema instances that don't conflict.
    //
    // Otherwise, create one DataAccess and hold on to the reference.
    if (this.databaseConfig.getBooleanOr("slick.createSchema", default = true)) {
      val schemaManager = SchemaManager.fromConfig(this.databaseConfig)
      val future = schemaManager.updateSchema(dataAccess.driver, dataAccess.schema, database)
      Await.result(future, Duration.Inf)
    }
  }

  /**
    * Create a special execution context, a fixed thread pool, to run each of our composite database actions. Running
    * each composite action as a runnable within the pool will ensure that-- at most-- the same number of actions are
    * running as there are available connections. Thus there should never be a connection deadlock, as outlined in
    * - https://github.com/slick/slick/issues/1274
    * - https://groups.google.com/d/msg/scalaquery/5MCUnwaJ7U0/NLLMotX9BQAJ
    *
    * Custom future thread pool based on:
    * - http://stackoverflow.com/questions/15285284/how-to-configure-a-fine-tuned-thread-pool-for-futures#comment23278672_15285441
    *
    * Database config parameter defaults based on: (expand the `forConfig` scaladoc for a full list of values)
    * - http://slick.typesafe.com/doc/3.1.0/api/index.html#slick.jdbc.JdbcBackend$DatabaseFactoryDef@forConfig(path:String,config:com.typesafe.config.Config,driver:java.sql.Driver,classLoader:ClassLoader):JdbcBackend.this.Database
    *
    * Reuses the error reporter from the database's executionContext.
    */
  private val actionThreadPool: ExecutorService = {
    val dbNumThreads = databaseConfig.getIntOr("db.numThreads", 20)
    val dbMaximumPoolSize = databaseConfig.getIntOr("db.maxConnections", dbNumThreads * 5)
    val actionThreadPoolSize = databaseConfig.getIntOr("actionThreadPoolSize", dbNumThreads) min dbMaximumPoolSize
    Executors.newFixedThreadPool(actionThreadPoolSize)
  }

  private val actionExecutionContext: ExecutionContext = ExecutionContext.fromExecutor(
    actionThreadPool, database.executor.executionContext.reportFailure
  )

  protected[this] lazy val useSlickUpserts =
    dataAccess.driver.capabilities.contains(JdbcProfile.capabilities.insertOrUpdate)

  protected[this] def assertUpdateCount(description: String, updates: Int, expected: Int): DBIO[Unit] = {
    if (updates == expected) {
      DBIO.successful(Unit)
    } else {
      DBIO.failed(new RuntimeException(s"$description expected update count $expected, got $updates"))
    }
  }

  override def close(): Unit = {
    actionThreadPool.shutdown()
    database.close()
  }

  protected[this] def runTransaction[R](action: DBIO[R]): Future[R] = {
    //database.run(action.transactionally) <-- https://github.com/slick/slick/issues/1274
    Future(Await.result(database.run(action.transactionally), Duration.Inf))(actionExecutionContext)
  }
}
