package cromwell.database.slick

import java.sql.{Connection, PreparedStatement, Statement}
import java.util.concurrent.{ExecutorService, Executors}

import com.mysql.cj.jdbc.exceptions.MySQLTransactionRollbackException
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.database.slick.tables.DataAccessComponent
import cromwell.database.sql.SqlDatabase
import net.ceedubs.ficus.Ficus._
import org.postgresql.util.{PSQLException, ServerErrorMessage}
import org.slf4j.LoggerFactory
import slick.basic.DatabaseConfig
import slick.jdbc.{JdbcCapabilities, JdbcProfile, PostgresProfile, TransactionIsolation}

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
    if (slickDatabase.databaseConfig.getOrElse("slick.createSchema", true)) {
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

  override lazy val connectionDescription = databaseConfig.getString(urlKey)

  SlickDatabase.log.info(s"Running with database $urlKey = $connectionDescription")

  /**
    * Original as of Slick 3.1.0:
    * ---------------------------
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
    *
    * Update as of Slick 3.2.3:
    * -------------------------
    * Even though https://github.com/slick/slick/issues/1274 has been marked resolved, we are still seeing sporadic
    * deadlocks in ServicesStoreSpec's `it should "not deadlock"` with Slick 3.2.3. The deadlocks all seem to appear
    * shortly after the logged error "count cannot be increased".
    *
    * Based on logs, this error does not appear in production, only in tests.
    * Also, "count cannot be decreased" is likely a **different** error.
    *
    * Among many possibilities:
    * - The test configuration is overly conservative with the two db threads and two db connections.
    * - Similarly, it's possible that by switching to this other thread pool for transactions it ends up freeing Slick's
    *   internal thread pool just enough to fix the problem.
    * - Or perhaps this doesn't fix the problem at all and we will see the failures in tests again soon.
    *
    * Further reading:
    * - https://github.com/broadinstitute/cromwell/issues/4328#issuecomment-439187592
    * - https://github.com/slick/slick/pull/1914
    * - https://github.com/slick/slick/issues/1614
    * - https://github.com/brettwooldridge/HikariCP/issues/1084
    * - https://github.com/slick/slick/issues/1807#issuecomment-358244898
    *
    * Database config parameter defaults updated to use the Slick 3.2.3 values:
    * (expand the `forConfig` scaladoc for a full list of values)
    * - http://slick.typesafe.com/doc/3.2.3/api/index.html#slick.jdbc.JdbcBackend$DatabaseFactoryDef@forConfig(path:String,config:com.typesafe.config.Config,driver:java.sql.Driver,classLoader:ClassLoader):JdbcBackend.this.Database
    *
    * Specifically, maxConnections has been lowered from (numThreads * 5) down to (numThreads * 1) to match:
    * https://github.com/slick/slick/commit/c79d50c55747515314f26b5f4749bb2b22084fe0
    */
  private val actionThreadPool: ExecutorService = {
    val dbNumThreads = databaseConfig.getOrElse("db.numThreads", 20)
    val dbMaximumPoolSize = databaseConfig.getOrElse("db.maxConnections", dbNumThreads)
    val actionThreadPoolSize = databaseConfig.getOrElse("actionThreadPoolSize", dbNumThreads) min dbMaximumPoolSize
    Executors.newFixedThreadPool(actionThreadPoolSize)
  }

  private val actionExecutionContext: ExecutionContext = ExecutionContext.fromExecutor(
    actionThreadPool, database.executor.executionContext.reportFailure
  )

  protected[this] lazy val insertBatchSize = databaseConfig.getOrElse("insert-batch-size", 2000)

  protected[this] lazy val useSlickUpserts =
    dataAccess.driver.capabilities.contains(JdbcCapabilities.insertOrUpdate)

  protected[this] def assertUpdateCount(description: String, updates: Int, expected: Int): DBIO[Unit] = {
    if (updates == expected) {
      DBIO.successful(())
    } else {
      DBIO.failed(new RuntimeException(s"$description expected update count $expected, got $updates"))
    }
  }

  override def withConnection[A](block: Connection => A): A = {
    /*
     TODO: Should this withConnection() method have a (implicit?) timeout parameter, that it passes on to Await.result?
     If we run completely asynchronously, nest calls to withConnection, and then call flatMap, the outer connection may
     already be closed before an inner block finishes running.
     */
    Await.result(database.run(SimpleDBIO(context => block(context.connection))), Duration.Inf)
  }

  override def close(): Unit = {
    actionThreadPool.shutdown()
    database.close()
  }

  private val debugExitConfigPath = "danger.debug.only.exit-on-rollback-exception-with-status-code"
  private val debugExitStatusCodeOption = ConfigFactory.load.getAs[Int](debugExitConfigPath)

  protected[this] def runTransaction[R](action: DBIO[R],
                                        isolationLevel: TransactionIsolation = TransactionIsolation.RepeatableRead,
                                        timeout: Duration = Duration.Inf): Future[R] = {
    runActionInternal(action.transactionally.withTransactionIsolation(isolationLevel), timeout = timeout)
  }

  /* Note that this is only appropriate for actions that do not involve Blob
   * or Clob fields in Postgres, since large object support requires running
   * transactionally.  Use runLobAction instead, which will still run in
   * auto-commit mode when using other database engines.
   */
  protected[this] def runAction[R](action: DBIO[R]): Future[R] = {
    runActionInternal(action.withPinnedSession)
  }

  /* Wrapper for queries where Clob/Blob types are used
   * https://stackoverflow.com/questions/3164072/large-objects-may-not-be-used-in-auto-commit-mode#answer-3164352
   */
  protected[this] def runLobAction[R](action: DBIO[R]): Future[R] = {
    dataAccess.driver match {
      case PostgresProfile => runTransaction(action)
      case _ => runAction(action)
    }
  }

  private def runActionInternal[R](action: DBIO[R], timeout: Duration = Duration.Inf): Future[R] = {
    //database.run(action) <-- See comment above private val actionThreadPool
    Future {
      try {
        if (timeout.isFinite()) {
          // https://stackoverflow.com/a/52569275/818054
          Await.result(database.run(action.withStatementParameters(statementInit = _.setQueryTimeout(timeout.toSeconds.toInt))), Duration.Inf)
        } else {
          Await.result(database.run(action), Duration.Inf)
        }
      } catch {
        case rollbackException: MySQLTransactionRollbackException =>
          debugExitStatusCodeOption match {
            case Some(status) =>
              SlickDatabase.log.error("Got a rollback!", rollbackException)
              System.err.println(s"Exiting with code $status as $debugExitConfigPath is set")
              System.exit(status)
            case _ => /* keep going */
          }
          throw rollbackException
        case pSQLException: PSQLException =>
          val detailOption = for {
            message <- Option(pSQLException.getServerErrorMessage)
            detail <- Option(message.getDetail)
          } yield detail

          detailOption match {
            case None => throw pSQLException
            case Some(_) =>
              /*
              The exception may contain possibly sensitive row contents within the DETAIL section. Remove it.

              Tried adjusting this using configuration:
              - log_error_verbosity=TERSE
              - log_min_messages=PANIC
              - client_min_messages=ERROR

              Instead resorting to reflection.
               */
              val message = pSQLException.getServerErrorMessage
              val field = classOf[ServerErrorMessage].getDeclaredField("m_mesgParts")
              field.setAccessible(true)
              val parts = field.get(message).asInstanceOf[java.util.Map[Character, String]]
              parts.remove('D')
              // The original exception has already stored the DETAIL into a string. So we must create a new Exception.
              throw new PSQLException(message)
          }
      }
    }(actionExecutionContext)
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
      context.session.withPreparedStatement[Array[Int]](compiled.upsert.sql) { st: PreparedStatement =>
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
