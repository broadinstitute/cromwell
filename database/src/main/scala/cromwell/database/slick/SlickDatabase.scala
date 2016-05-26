package cromwell.database.slick

import java.sql.{Clob, Timestamp}
import java.util.UUID
import java.util.concurrent.{ExecutorService, Executors}

import com.typesafe.config.{Config, ConfigFactory, ConfigValueFactory}
import cromwell.database.SqlDatabase
import cromwell.database.obj._
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
      * "${slick.uniqueSchema}".
      *
      * This allows each instance of a SlickDataAccess object to use a clean, and different, in memory database.
      *
      * @return Config with ${slick.uniqueSchema} in url replaced with a unique string.
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
class SlickDatabase(databaseConfig: Config) extends SqlDatabase {

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

  private lazy val useSlickUpserts = dataAccess.driver.capabilities.contains(JdbcProfile.capabilities.insertOrUpdate)

  override def close(): Unit = {
    actionThreadPool.shutdown()
    database.close()
  }

  private def runTransaction[R](action: DBIO[R]): Future[R] = {
    //database.run(action.transactionally) <-- https://github.com/slick/slick/issues/1274
    Future(Await.result(database.run(action.transactionally), Duration.Inf))(actionExecutionContext)
  }

  override def createWorkflow(workflowExecution: WorkflowExecution,
                              workflowExecutionAuxFromWorkflowExecutionId: Int => WorkflowExecutionAux,
                              workflowInputsFromWorkflowExecutionId: Int => Seq[Symbol],
                              workflowExecutionsFromWorkflowExecutionId: Int => Seq[Execution],
                              workflowExecutionInfosFromExecutionId: Int => Seq[ExecutionInfo])
                             (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      workflowExecutionId <- dataAccess.workflowExecutionIdsAutoInc += workflowExecution
      _ <- dataAccess.workflowExecutionAuxIdsAutoInc += workflowExecutionAuxFromWorkflowExecutionId(workflowExecutionId)
      _ <- dataAccess.symbolIdsAutoInc ++= workflowInputsFromWorkflowExecutionId(workflowExecutionId)
      // Insert an execution row
      executionIds <- dataAccess.executionIdsAutoInc ++= workflowExecutionsFromWorkflowExecutionId(workflowExecutionId)
      _ <- dataAccess.executionInfoIdsAutoInc ++= executionIds flatMap workflowExecutionInfosFromExecutionId
    } yield ()

    runTransaction(action)
  }

  override def getWorkflowExecution(workflowUuid: String)
                                   (implicit ec: ExecutionContext): Future[WorkflowExecution] = {
    val action = dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowUuid).result.head

    runTransaction(action)
  }

  override def upsertRuntimeAttributes(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                       attributes: Map[String, String])
                                      (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      executionId <- dataAccess.executionIdsByWorkflowExecutionUuidAndCallKey(
        workflowUuid, callFqn, index, attempt).result.head
      _ <- DBIO.sequence(attributes map upsertRuntimeAttribute(executionId))
    } yield ()

    runTransaction(action)
  }

  private def upsertRuntimeAttribute(executionId: Int)(keyValue: (String, String))
                                    (implicit ec: ExecutionContext): DBIO[Unit] = {
    val (key, value) = keyValue

    if (useSlickUpserts) {
      for {
        _ <- dataAccess.runtimeAttributeIdsAutoInc.insertOrUpdate(RuntimeAttribute(executionId, key, value))
      } yield ()
    } else {
      for {
        updateCount <- dataAccess.runtimeAttributeValuesByExecutionIdAndName(executionId, key).update(value)
        _ <- updateCount match {
          case 0 => dataAccess.runtimeAttributeIdsAutoInc += RuntimeAttribute(executionId, key, value)
          case 1 => DBIO.successful(Unit)
          case _ => DBIO.failed(new RuntimeException(s"Unexpected backend update count $updateCount"))
        }
      } yield ()
    }
  }

  override def getAllRuntimeAttributes(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[Traversable[(String, Int, Int, String, String)]] = {
    val action = dataAccess.runtimeAttributesByWorkflowUUID(workflowUuid).result

    runTransaction(action)
  }

  override def insertCalls(workflowUuid: String,
                           executionsFromWorkflowExecutionId: Int => Seq[Execution],
                           executionInfosFromWorkflowExecutionId: Int => Seq[ExecutionInfo])
                          (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      workflowExecutionId <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      executionIds <- dataAccess.executionIdsAutoInc ++= executionsFromWorkflowExecutionId(workflowExecutionId)
      // Add the empty execution info rows using the execution as the FK.
      _ <- dataAccess.executionInfoIdsAutoInc ++= executionIds flatMap executionInfosFromWorkflowExecutionId
    } yield ()

    runTransaction(action)
  }

  override def getWorkflowStateString(workflowUuid: String)
                                     (implicit ec: ExecutionContext): Future[Option[String]] = {
    val action = dataAccess.workflowExecutionStatusByWorkflowExecutionUuid(workflowUuid).result.headOption

    runTransaction(action)
  }

  override def getExecutionStatuses(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[Traversable[(String, Int, Int, String, Option[Int], Option[String], Option[String])]] = {
    val action = for {
    // NOTE: For now, intentionally causes query to error out instead of returning an Map.empty
      _ <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      result <- dataAccess.executionCallStatusesByWorkflowExecutionUuid(workflowUuid).result
    } yield result

    runTransaction(action)
  }

  override def getExecutionStatuses(workflowUuid: String, callFqn: String)(implicit ec: ExecutionContext):
  Future[Traversable[(String, Int, Int, String, Option[Int], Option[String], Option[String])]] = {
    val action = dataAccess.executionCallStatusesByWorkflowExecutionUuidAndCallFqn(workflowUuid, callFqn).result

    runTransaction(action)
  }

  override def getExecutionStatus(workflowUuid: String, callFqn: String, index: Int, attempt: Int)
                                 (implicit ec: ExecutionContext):
  Future[Option[(String, Option[Int], Option[String], Option[String])]] = {
    val action = dataAccess.executionCallStatusesByWorkflowExecutionUuidAndCallKey(
      workflowUuid, callFqn, index, attempt).result.headOption

    runTransaction(action)
  }

  override def getWorkflowExecutionAndAuxTuple(workflowExecutionId: Int)(implicit ec: ExecutionContext):
  Future[(WorkflowExecution, WorkflowExecutionAux)] = {
    val action = dataAccess.workflowExecutionsAndAuxesByWorkflowExecutionId(workflowExecutionId).result.head

    runTransaction(action)
  }

  override def getWorkflowExecutionAndAuxTuple(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[(WorkflowExecution, WorkflowExecutionAux)] = {
    val action = dataAccess.workflowExecutionsAndAuxesByWorkflowExecutionUuid(workflowUuid).result.head

    runTransaction(action)
  }

  override def getWorkflowExecutionAndAuxTuples(states: Traversable[String])(implicit ec: ExecutionContext):
  Future[Traversable[(WorkflowExecution, WorkflowExecutionAux)]] = {
    val action = dataAccess.workflowExecutionAndAuxesByStatuses(states).result

    runTransaction(action)
  }

  override def updateWorkflowState(workflowUuid: String, workflowState: String, endDate: Option[Timestamp])
                                  (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      updateCount <- dataAccess.
        workflowExecutionStatusEndDtByWorkflowExecutionUuid(workflowUuid).
        update(workflowState, endDate)
      _ <- assertUpdateCount("updateWorkflowState", updateCount, 1)
    } yield ()

    runTransaction(action)
  }

  override def getExecutionInfos(workflowUuid: String, callFqn: String, attempt: Int)
                                (implicit ec: ExecutionContext): Future[Traversable[ExecutionInfo]] = {
    val action = dataAccess.
      executionInfosByWorkflowExecutionUuidAndCallFqnAndAttempt(workflowUuid, callFqn, attempt).result

    runTransaction(action)
  }

  override def getExecutionInfoByKey(workflowUuid: String, callFqn: String, attempt: Int, key: String)
                                    (implicit ec: ExecutionContext): Future[Option[Option[String]]] = {
    val action = dataAccess.executionInfoValueByWorkflowExecutionUuidAndCallFqnAndAttemptAndKey(
      workflowUuid, callFqn, attempt, key).result.headOption

    runTransaction(action)
  }

  override def updateExecutionInfo(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                   key: String, value: Option[String])(implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      executionId <- dataAccess.executionIdsByWorkflowExecutionUuidAndCallKey(
        workflowUuid, callFqn, index, attempt).result.head

      updateCount <- dataAccess.executionInfoValueByExecutionAndKey(executionId, key).update(value)
      _ <- assertUpdateCount("updateExecutionInfo", updateCount, 1)
    } yield ()

    runTransaction(action)
  }

  override def upsertExecutionInfo(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                   keyValues: Map[String, Option[String]])
                                  (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      executionId <- dataAccess.executionIdsByWorkflowExecutionUuidAndCallKey(
        workflowUuid, callFqn, index, attempt).result.head
      _ <- DBIO.sequence(keyValues map upsertExecutionInfo(executionId))
    } yield ()

    runTransaction(action)
  }

  protected def upsertExecutionInfo(executionId: Int)(keyValue: (String, Option[String]))
                                   (implicit ec: ExecutionContext): DBIO[Unit] = {
    val (key, value) = keyValue

    if (useSlickUpserts) {
      for {
        _ <- dataAccess.executionInfoIdsAutoInc.insertOrUpdate(ExecutionInfo(executionId, key, value))
      } yield ()
    } else {
      for {
        count <- dataAccess.executionInfoValueByExecutionAndKey(executionId, key).update(value)
        _ <- count match {
          case 0 => dataAccess.executionInfoIdsAutoInc += ExecutionInfo(executionId, key, value)
          case 1 => DBIO.successful(Unit)
          case _ => DBIO.failed(new RuntimeException(s"Unexpected update count $count"))
        }
      } yield ()
    }
  }

  override def getAllSymbols(workflowUuid: String)
                            (implicit ec: ExecutionContext): Future[Traversable[Symbol]] = {
    val action = dataAccess.symbolsByWorkflowExecutionUuid(workflowUuid).result

    runTransaction(action)
  }

  /** Returns all NON SHARDS outputs for this workflowId */
  override def getWorkflowOutputSymbols(workflowUuid: String)
                                       (implicit ec: ExecutionContext): Future[Traversable[Symbol]] = {
    val action = dataAccess.symbolsForWorkflowOutput(workflowUuid).result

    runTransaction(action)
  }

  override def getAllSymbols(workflowUuid: String, ioValue: String)
                            (implicit ec: ExecutionContext): Future[Traversable[Symbol]] = {
    val action = dataAccess.symbolsByWorkflowExecutionUuidAndIo(workflowUuid, ioValue).result

    runTransaction(action)
  }

  override def getAllSymbols(workflowUuid: String, ioValue: String, callFqn: String, index: Int)
                            (implicit ec: ExecutionContext): Future[Traversable[Symbol]] = {
    val action = dataAccess.symbolsByWorkflowExecutionUuidAndIoAndScopeAndIndex(
      workflowUuid, ioValue, callFqn, index).result

    runTransaction(action)
  }

  override def setOutputs(workflowUuid: String, symbolsFromWorkflowExecutionId: Int => Seq[Symbol])
                         (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      workflowExecutionId <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      _ <- dataAccess.symbolIdsAutoInc ++= symbolsFromWorkflowExecutionId(workflowExecutionId)
    } yield ()

    runTransaction(action)
  }

  /**
    * Updates the existing input symbols to replace expressions with real values.
    *
    * @return The number of rows updated per callInput - as a Future.
    */
  override def updateCallInputs(workflowUuid: String, callFqn: String, index: Int,
                                callInputs: Traversable[(String, String, Option[Clob])])
                               (implicit ec: ExecutionContext): Future[Traversable[Int]] = {
    val action = for {
      workflowExecutionId <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      updates = callInputs map { case (inputName, wdlType, wdlValue) =>
        dataAccess.symbolWdlTypeAndWdlValueByWorkflowAndScopeAndIndexAndName(
          workflowExecutionId, callFqn, index, inputName).update(wdlType, wdlValue)
      }
      count <- DBIO.sequence(updates)
    } yield count

    runTransaction(action)
  }

  override def setExecutionEvents(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                  executionEventsFromExecutionId: Int => Seq[ExecutionEvent])
                                 (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      executionId <- dataAccess.executionIdsByWorkflowExecutionUuidAndCallKey(
        workflowUuid, callFqn, index, attempt).result.head
      _ <- dataAccess.executionEventIdsAutoInc ++= executionEventsFromExecutionId(executionId)
    } yield ()

    runTransaction(action)
  }

  override def setExecutionEvents(workflowUuid: String, callFqn: String, attempt: Int,
                                  executionEventsFromExecutionId: Int => Seq[ExecutionEvent])
                                 (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      executionId <- dataAccess.executionIdsByWorkflowExecutionUuidAndCallFqnAndAttempt(
        workflowUuid, callFqn, attempt).result.head
      _ <- dataAccess.executionEventIdsAutoInc ++= executionEventsFromExecutionId(executionId)
    } yield ()

    runTransaction(action)
  }

  override def getAllExecutionEvents(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[Traversable[(String, Int, Int, String, Timestamp, Timestamp)]] = {
    val action = dataAccess.executionEventsByWorkflowExecutionUuid(workflowUuid).result

    runTransaction(action)
  }

  override def addCallFailureEvent(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                   failureEventsFromWorkflowExecutionIdAndExecutionId: (Int, Int) => FailureEvent)
                                  (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      (executionId, workflowExecutionId) <- dataAccess.
        executionIdAndWorkflowExecutionIdsByWorkflowExecutionUuidAndCallKey(
          workflowUuid, callFqn, index, attempt).result.head
      _ <- dataAccess.failureEventIdsAutoInc +=
        failureEventsFromWorkflowExecutionIdAndExecutionId(workflowExecutionId, executionId)
    } yield ()

    runTransaction(action)
  }

  override def addCallFailureEvent(workflowUuid: String, callFqn: String, attempt: Int,
                                   failureEventsFromWorkflowExecutionIdAndExecutionId: (Int, Int) => FailureEvent)
                                  (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      (executionId, workflowExecutionId) <- dataAccess.
        executionIdAndWorkflowExecutionIdsByWorkflowExecutionUuidAndCallFqnAndAttempt(
          workflowUuid, callFqn, attempt).result.head
      _ <- dataAccess.failureEventIdsAutoInc +=
        failureEventsFromWorkflowExecutionIdAndExecutionId(workflowExecutionId, executionId)
    } yield ()

    runTransaction(action)
  }

  override def addWorkflowFailureEvent(workflowUuid: String, failureEventsFromWorkflowExecutionId: Int => FailureEvent)
                                      (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      workflowExecutionId <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      _ <- dataAccess.failureEventIdsAutoInc += failureEventsFromWorkflowExecutionId(workflowExecutionId)
    } yield ()

    runTransaction(action)
  }

  override def getFailureEvents(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[Traversable[(String, String, Timestamp, Option[String], Option[Int], Option[Int])]] = {
    val action = dataAccess.failureEventsByWorkflowExecutionUuid(workflowUuid).result

    runTransaction(action)
  }

  private def assertUpdateCount(method: String, updates: Int, expected: Int): DBIO[Unit] = {
    if (updates == expected) {
      DBIO.successful(Unit)
    } else {
      DBIO.failed(new RuntimeException(s"$method expected update count $expected, got $updates"))
    }
  }

  override def setStartingStatus(workflowUuid: String, statusString: String, startDt: Option[Timestamp],
                                 scopeKeys: Traversable[(String, Int, Int)])
                                (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      workflowExecutionId <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      updates = scopeKeys map {
        case (callFqn, index, attempt) => for {
          updateCount <-
          dataAccess.
            executionStatusesAndStartDtByWorkflowExecutionIdAndCallKey(workflowExecutionId, callFqn, index, attempt).
            update(statusString, startDt)
          _ = assertUpdateCount("updateStatus", updateCount, 1)
        } yield ()
      }
      _ <- DBIO.sequence(updates)
    } yield ()

    runTransaction(action)
  }

  override def updateStatus(workflowUuid: String, statusString: String, scopeKeys: Traversable[(String, Int, Int)])
                           (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      workflowExecutionId <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      updates = scopeKeys map {
        case (callFqn, index, attempt) => for {
          updateCount <- dataAccess.
            executionStatusesByWorkflowExecutionIdAndCallKey(workflowExecutionId, callFqn, index, attempt).
            update(statusString)
          _ = assertUpdateCount("updateStatus", updateCount, 1)
        } yield ()
      }
      _ <- DBIO.sequence(updates)
    } yield ()

    runTransaction(action)
  }

  override def setTerminalWithoutClone(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                       statusString: String, endDt: Option[Timestamp], scriptReturnCode: Option[Int],
                                       overallHash: Option[String], dockerHash: Option[String])
                                      (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      workflowExecutionId <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      updateCount <- dataAccess.
        executionsByWorkflowExecutionIdAndCallKey(workflowExecutionId, callFqn, index, attempt).
        update(statusString, endDt, scriptReturnCode, overallHash, dockerHash, None)
      _ <- assertUpdateCount("setTerminalWithoutClone", updateCount, 1)
    } yield ()

    runTransaction(action)
  }

  override def setTerminalWithClone(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                    statusString: String, endDt: Option[Timestamp], scriptReturnCode: Option[Int],
                                    overallHash: Option[String], dockerHash: Option[String],
                                    workflowUuidClone: String, callFqnClone: String, indexClone: Int,
                                    attemptClone: Int)(implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      clonedFromId <- dataAccess.executionIdsByWorkflowExecutionUuidAndCallKey(
        workflowUuidClone, callFqnClone, indexClone, attemptClone).result.headOption
      _ <- clonedFromId match {
        case Some(_) => DBIO.successful(Unit)
        case None => DBIO.failed(new RuntimeException(
          s"Clone from not found for ($workflowUuidClone, $callFqnClone, $indexClone, $attemptClone)"))
      }
      workflowExecutionId <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      updateCount <- dataAccess.
        executionsByWorkflowExecutionIdAndCallKey(workflowExecutionId, callFqn, index, attempt).
        update(statusString, endDt, scriptReturnCode, overallHash, dockerHash, clonedFromId)
      _ <- assertUpdateCount("setTerminalWithClone", updateCount, 1)
    } yield ()

    runTransaction(action)
  }

  override def getExecutions(workflowUuid: String)
                            (implicit ec: ExecutionContext): Future[Traversable[Execution]] = {
    val action = dataAccess.executionsByWorkflowExecutionUuid(workflowUuid).result

    runTransaction(action)
  }

  override def getExecutionsWithResuableResultsByHash(hash: String)
                                                     (implicit ec: ExecutionContext): Future[Traversable[Execution]] = {
    val action = dataAccess.executionsWithReusableResultsByExecutionHash(hash).result

    runTransaction(action)
  }

  override def updateWorkflowOptions(workflowUuid: String, workflowOptionsJson: Clob)
                                    (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      workflowExecutionId <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      updateCount <- dataAccess.workflowOptionsByWorkflowExecutionId(workflowExecutionId).update(workflowOptionsJson)
      _ <- assertUpdateCount("updateWorkflowOptions", updateCount, 1)
    } yield ()

    runTransaction(action)
  }

  override def queryWorkflowExecutions(statuses: Set[String], names: Set[String], uuids: Set[String],
                                       startDate: Option[Timestamp], endDate: Option[Timestamp],
                                       page: Option[Int], pageSize: Option[Int])
                                      (implicit ec: ExecutionContext): Future[Traversable[WorkflowExecution]] = {
    val action = dataAccess.queryWorkflowExecutions(statuses, names, uuids, startDate, endDate, page, pageSize).result

    runTransaction(action)
  }

  override def countWorkflowExecutions(statuses: Set[String], names: Set[String], uuids: Set[String],
                                       startDate: Option[Timestamp], endDate: Option[Timestamp])
                                      (implicit ec: ExecutionContext): Future[Int] = {
    val action = dataAccess.countWorkflowExecutions(statuses, names, uuids, startDate, endDate).result

    runTransaction(action)
  }

  override def updateCallCaching(workflowUuid: String, allow: Boolean)(implicit ec: ExecutionContext): Future[Int] = {
    val action = for {
      workflowExecutionId <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      count <- dataAccess.executionAllowsResultReusesByWorkflowExecutionId(workflowExecutionId).update(allow)
    } yield count

    runTransaction(action)
  }

  override def updateCallCaching(workflowUuid: String, allow: Boolean, callFqn: String, attempt: Int)
                                (implicit ec: ExecutionContext): Future[Int] = {
    val action = for {
      workflowExecutionId <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      count <- dataAccess.
        executionAllowsResultReusesByWorkflowExecutionIdAndCallFqnAndAttempt(workflowExecutionId, callFqn, attempt).
        update(allow)
    } yield count

    runTransaction(action)
  }

  override def updateCallCaching(workflowUuid: String, allow: Boolean, callFqn: String, index: Int, attempt: Int)
                                (implicit ec: ExecutionContext): Future[Int] = {
    val action = for {
      workflowExecutionId <- dataAccess.workflowExecutionIdsByWorkflowExecutionUuid(workflowUuid).result.head
      count <- dataAccess.executionAllowsResultReusesByWorkflowExecutionIdAndCallKey(
        workflowExecutionId, callFqn, index, attempt).update(allow)
    } yield count

    runTransaction(action)
  }

  override def infosByExecution(workflowUuid: String)
                               (implicit ec: ExecutionContext): Future[Traversable[(Execution, ExecutionInfo)]] = {
    val action = dataAccess.executionsAndExecutionInfosByWorkflowExecutionUuid(workflowUuid).result

    runTransaction(action)
  }

  override def infosByExecution(workflowUuid: String, callFqn: String)
                               (implicit ec: ExecutionContext): Future[Traversable[(Execution, ExecutionInfo)]] = {
    val action = dataAccess.executionsAndExecutionInfosByWorkflowExecutionUuidAndCallFqn(workflowUuid, callFqn).result

    runTransaction(action)
  }

  override def callCacheDataByExecution(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[Traversable[(Execution, Option[String], Option[String])]] = {
    val action = dataAccess.selectExecutionsWithCacheHitWorkflowAndCall(workflowUuid)

    runTransaction(action)
  }

  override def runningExecutionsAndExecutionInfos(workflowUuid: String, statuses: Set[String])(implicit ec: ExecutionContext):
  Future[Traversable[(Execution, ExecutionInfo)]] = {
    val action = dataAccess.runningExecutionsAndExecutionInfosByWorkflowExecutionUuid(workflowUuid, statuses).result

    runTransaction(action)
  }

  override def addMetadataEvent(workflowUuid: String,
                                key: String,
                                value: String,
                                valueType: String,
                                timestamp: Timestamp)(implicit ec: ExecutionContext): Future[Unit] = {

    val action = dataAccess.metadataAutoInc += Metadatum(workflowUuid, key,
      callFqn = None, index = None, attempt = None, Option(value), Option(valueType), timestamp)
    runTransaction(action).map(_ => ())
  }

  override def addMetadataEvent(workflowUuid: String,
                                key: String,
                                callFqn: String,
                                index: Option[Int],
                                attempt: Int,
                                value: String,
                                valueType: String,
                                timestamp: Timestamp)(implicit ec: ExecutionContext): Future[Unit] = {

    val action = dataAccess.metadataAutoInc += Metadatum(workflowUuid, key,
      Option(callFqn), index, Option(attempt), Option(value), Option(valueType), timestamp)
    runTransaction(action).map(_ => ())
  }

  override def queryMetadataEvents(workflowUuid: String)
                                  (implicit ec: ExecutionContext): Future[Seq[Metadatum]] = {

    val action = dataAccess.metadataByWorkflowUuid(workflowUuid).result
    runTransaction(action)
  }

  override def queryMetadataEvents(workflowUuid: String,
                                   key: String)
                                  (implicit ec: ExecutionContext): Future[Seq[Metadatum]] = {

    val action = dataAccess.metadataByWorkflowUuidAndKey(workflowUuid, key).result
    runTransaction(action)
  }

  override def queryMetadataEvents(workflowUuid: String,
                                   callFqn: String,
                                   index: Option[Int],
                                   attempt: Int)
                                  (implicit ec: ExecutionContext): Future[Seq[Metadatum]] = {

    val action = dataAccess.metadataByWorkflowUuidAndCallFqnAndIndexAndAttempt(workflowUuid, callFqn, index, attempt).result
    runTransaction(action)
  }

  override def queryMetadataEvents(workflowUuid: String,
                                   key: String,
                                   callFqn: String,
                                   index: Option[Int],
                                   attempt: Int)
                                  (implicit ec: ExecutionContext): Future[Seq[Metadatum]] = {

    val action = dataAccess.metadataByWorkflowUuidAndKeyAndCallFqnAndIndexAndAttempt(workflowUuid, key, callFqn, index, attempt).result
    runTransaction(action)
  }
}
