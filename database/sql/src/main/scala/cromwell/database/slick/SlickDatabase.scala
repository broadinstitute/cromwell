package cromwell.database.slick

import java.sql.{Connection, PreparedStatement, Statement}

import com.typesafe.config.Config
import cromwell.database.slick.tables.DataAccessComponent
import cromwell.database.sql.SqlDatabase
import net.ceedubs.ficus.Ficus._
import org.slf4j.LoggerFactory
import slick.basic.DatabaseConfig
import slick.jdbc.{JdbcCapabilities, JdbcProfile}

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, Future}

object SlickDatabase {
  /**
    * Returns either the "url" or "properties.url"
    */
  def urlKey(config: Config) = if (config.hasPath("db.url")) "db.url" else "db.properties.url"

  lazy val log = LoggerFactory.getLogger("cromwell.database.slick")

  def createSchema(slickDatabase: SlickDatabase): Unit = {
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
    // The value `${uniqueSchema}` may be used in the url, in combination with `slick.createSchema = true`, to
    // generate unique schema instances that don't conflict.
    //
    // Otherwise, create one DataAccess and hold on to the reference.
    if (slickDatabase.databaseConfig.as[Option[Boolean]]("slick.createSchema").getOrElse(true)) {
      import slickDatabase.dataAccess.driver.api._
      Await.result(slickDatabase.database.run(slickDatabase.dataAccess.schema.create), Duration.Inf)
    }
  }

  def getDatabaseConfig(name: String, parentConfig: Config): Config = {
    val rootDatabaseConfig = parentConfig.getConfig("database")
    rootDatabaseConfig.getOrElse(name, rootDatabaseConfig)
  }
}

/**
  * Data Access implementation using Slick.
  */
abstract class SlickDatabase(override val originalDatabaseConfig: Config) extends SqlDatabase {

  override val urlKey = SlickDatabase.urlKey(originalDatabaseConfig)
  protected val slickConfig = DatabaseConfig.forConfig[JdbcProfile]("", databaseConfig)

  /*
  Not a def because we need to have a "stable identifier" for the imports below.
  Must be overridden as a lazy val, or "early definitions", to avoid getting nulls on init.
  http://docs.scala-lang.org/tutorials/FAQ/initialization-order.html
   */
  val dataAccess: DataAccessComponent

  // Allows creation of a Database, plus implicits for running transactions
  import dataAccess.driver.api._

  // NOTE: if you want to refactor database is inner-class type: this.dataAccess.driver.backend.DatabaseFactory
  val database = slickConfig.db

  SlickDatabase.log.info(s"Running with database $urlKey = ${databaseConfig.getString(urlKey)}")

  protected[this] lazy val insertBatchSize = databaseConfig.as[Option[Int]]("insert-batch-size").getOrElse(2000)

  protected[this] lazy val useSlickUpserts =
    dataAccess.driver.capabilities.contains(JdbcCapabilities.insertOrUpdate)

  protected[this] def assertUpdateCount(description: String, updates: Int, expected: Int): DBIO[Unit] = {
    if (updates == expected) {
      DBIO.successful(())
    } else {
      DBIO.failed(new RuntimeException(s"$description expected update count $expected, got $updates"))
    }
  }

  override def withConnection[A](block: (Connection) => A): A = {
    /*
     TODO: Should this withConnection() method have a (implicit?) timeout parameter, that it passes on to Await.result?
     If we run completely asynchronously, nest calls to withConnection, and then call flatMap, the outer connection may
     already be closed before an inner block finishes running.
     */
    Await.result(database.run(SimpleDBIO(context => block(context.connection))), Duration.Inf)
  }

  override def close(): Unit = {
    database.close()
  }

  protected[this] def runTransaction[R](action: DBIO[R]): Future[R] = {
    database.run(action.transactionally)
  }

  /*
    * Upsert the provided values in batch.
    * Fails the query if one or more upsert failed.
    * Adapted from https://github.com/slick/slick/issues/1781
   */
  protected[this] def createBatchUpsert[T](description: String,
                                           compiled: dataAccess.driver.JdbcCompiledInsert,
                                           values: Iterable[T]
                                          )(implicit ec: ExecutionContext): DBIO[Unit] = {
    SimpleDBIO { context =>
      context.session.withPreparedStatement[Array[Int]](compiled.upsert.sql) { (st: PreparedStatement) =>
        values.foreach { update =>
          st.clearParameters()
          compiled.upsert.converter.set(update, st)
          st.addBatch()
        }
        st.executeBatch()
      }
    } flatMap { upsertCounts =>
      val failures = upsertCounts.filter(_ == Statement.EXECUTE_FAILED)
      if (failures.isEmpty) DBIO.successful(())
      else {
        val valueList = values.toList
        val failedRequests = failures.map(valueList(_))
        DBIO.failed(new RuntimeException(
          s"$description failed to upsert the following rows: ${failedRequests.mkString(", ")}"
        ))
      }
    }
  }
}
