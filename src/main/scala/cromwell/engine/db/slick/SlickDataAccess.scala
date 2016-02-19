package cromwell.engine.db.slick

import java.sql.{Clob, Timestamp}
import java.util.concurrent.{ExecutorService, Executors}
import java.util.{Date, UUID}
import javax.sql.rowset.serial.SerialClob

import _root_.slick.backend.DatabaseConfig
import _root_.slick.driver.JdbcProfile
import com.typesafe.config.{Config, ConfigFactory, ConfigValueFactory}
import cromwell.engine.ExecutionIndex._
import cromwell.engine.ExecutionStatus._
import cromwell.engine.backend.{BackendCall, Backend, JobKey, WorkflowQueryResult}
import cromwell.engine.db.DataAccess.ExecutionKeyToJobKey
import cromwell.engine.db._
import cromwell.engine.finalcall.FinalCall
import cromwell.engine.workflow._
import cromwell.engine.{CallOutput, SymbolHash, WorkflowOutputs, _}
import cromwell.webservice.{CallCachingParameters, WorkflowQueryParameters, WorkflowQueryResponse}
import lenthall.config.ScalaConfig._
import org.joda.time.DateTime
import org.slf4j.LoggerFactory
import wdl4s._
import wdl4s.types.{WdlPrimitiveType, WdlType}
import wdl4s.values.WdlValue

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, ExecutionContext, Future}
import scala.language.{implicitConversions, postfixOps}

object SlickDataAccess {
  type IoValue = String
  val IoInput = "INPUT"
  val IoOutput = "OUTPUT"

  lazy val rootConfig = ConfigFactory.load()
  private lazy val rootDatabaseConfig = rootConfig.getConfig("database")
  private lazy val databaseConfigName = rootDatabaseConfig.getStringOption("config")
  lazy val defaultDatabaseConfig = databaseConfigName.map(getDatabaseConfig).getOrElse(rootDatabaseConfig)

  def getDatabaseConfig(path: String) = rootDatabaseConfig.getConfig(path)

  implicit class DateToTimestamp(val date: Date) extends AnyVal {
    def toTimestamp = new Timestamp(date.getTime)
  }

  implicit class ClobToRawString(val clob: Clob) extends AnyVal {
    def toRawString: String = clob.getSubString(1, clob.length.toInt) // yes, it starts at 1
  }

  implicit class StringToClob(val str: String) extends AnyVal {
    def toClob: Clob = new SerialClob(str.toCharArray)
  }

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

  lazy val log = LoggerFactory.getLogger("slick")
}

/**
 * Data Access implementation using Slick.
 *
 * NOTE: the uses of .head below will cause an exception to be thrown
 * if the list is empty.  In every use case as of the writing of this comment,
 * those exceptions would have been wrapped in a failed Future and returned.
 */
class SlickDataAccess(databaseConfig: Config) extends DataAccess {

  import SlickDataAccess._

  def this() = this(SlickDataAccess.defaultDatabaseConfig)

  private val configWithUniqueSchema = this.databaseConfig.withUniqueSchema

  val slickConfig = DatabaseConfig.forConfig[JdbcProfile]("", configWithUniqueSchema)
  val dataAccess = new DataAccessComponent(slickConfig.driver)

  // Allows creation of a Database, plus implicits for running transactions
  import dataAccess.driver.api._

  // NOTE: if you want to refactor database is inner-class type: this.dataAccess.driver.backend.DatabaseFactory
  val database = slickConfig.db

  // Possibly create the database
  {
    import SlickDataAccess._
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

  private def wdlValueToDbValue(v: WdlValue): String = v.wdlType match {
    case p: WdlPrimitiveType => v.valueString
    case o => v.toWdlString
  }

  private def dbEntryToWdlValue(dbValue: String, wdlType: WdlType): WdlValue = wdlType match {
    // .get here is because we trust the value in the database is coercible to the given type
    case p: WdlPrimitiveType => p.coerceRawValue(dbValue).get
    case o => wdlType.fromWdlString(dbValue)
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
    actionThreadPool, database.executor.executionContext.reportFailure)

  override def close(): Unit = {
    actionThreadPool.shutdown()
    database.close()
  }

  // Run action with an outer transaction
  private def runTransaction[R](action: DBIOAction[R, _ <: NoStream, _ <: Effect]): Future[R] = {
    //database.run(action.transactionally) <-- https://github.com/slick/slick/issues/1274
    Future(Await.result(database.run(action.transactionally), Duration.Inf))(actionExecutionContext)
  }

  /**
   * Creates a row in each of the backend-info specific tables for each key in `keys` corresponding to the backend
   * `backend`.  Or perhaps defer this?
   */
  override def createWorkflow(workflowDescriptor: WorkflowDescriptor,
                              workflowInputs: Traversable[SymbolStoreEntry],
                              scopes: Traversable[Scope],
                              backend: Backend)(implicit ec: ExecutionContext): Future[Unit] = {

    val scopeKeys: Traversable[ExecutionStoreKey] = scopes collect {
      case call: Call => BackendCallKey(call, None, 1)
      case scatter: Scatter => ScatterKey(scatter, None)
      case finalCall: FinalCall => FinalCallKey(finalCall)
    }

    val action = for {

      workflowExecutionInsert <- dataAccess.workflowExecutionsAutoInc +=
        new WorkflowExecution(
          workflowDescriptor.id.toString,
          workflowDescriptor.name,
          WorkflowSubmitted.toString,
          new Date().toTimestamp)

      _ <- dataAccess.workflowExecutionAuxesAutoInc += new WorkflowExecutionAux(
        workflowExecutionInsert.workflowExecutionId.get,
        workflowDescriptor.sourceFiles.wdlSource.toClob,
        workflowDescriptor.sourceFiles.inputsJson.toClob,
        workflowDescriptor.sourceFiles.workflowOptionsJson.toClob
      )

      symbolInsert <- dataAccess.symbolsAutoInc ++= toInputSymbols(workflowExecutionInsert, workflowDescriptor.namespace.workflow, workflowInputs)

      // NOTE: Don't use DBIO.seq for **transforming** sequences
      // - DBIO.seq(mySeq: _*) runs *any* items in sequence, but converts Seq[ DBIOAction[_] ] to DBIOAction[ Unit ]
      // - DBIO.sequence(mySeq) converts Seq[ DBIOAction[R] ] to DBIOAction[ Seq[R] ]
      // - DBIO.fold(mySeq, init) converts Seq[ DBIOAction[R] ] to DBIOAction[R]

      _ <- DBIO.sequence(toScopeActions(workflowExecutionInsert, backend, scopeKeys))

    } yield ()

    runTransaction(action)
  }

  // Converts the SymbolStoreEntry to Symbols. Does not create the action to do the insert.
  private def toInputSymbols(workflowExecution: WorkflowExecution,
                             rootWorkflowScope: Workflow,
                             symbolStoreEntries: Traversable[SymbolStoreEntry]): Seq[Symbol] = {
    symbolStoreEntries.toSeq map { symbol =>
      val reportableResult = rootWorkflowScope.outputs exists { _.fullyQualifiedName == symbol.key.fqn }
      new Symbol(
        workflowExecution.workflowExecutionId.get,
        symbol.key.scope,
        symbol.key.name,
        symbol.key.index.fromIndex,
        if (symbol.key.input) IoInput else IoOutput,
        reportableResult,
        symbol.wdlType.toWdlString,
        symbol.wdlValue.map(v => wdlValueToDbValue(v).toClob),
        symbol.symbolHash map { _.value }
      )
    }
  }

  // Converts the Traversable[Call] to Seq[DBIOAction[]] that insert the correct rows
  private def toScopeActions(workflowExecution: WorkflowExecution, backend: Backend,
                            keys: Traversable[ExecutionStoreKey])(implicit ec: ExecutionContext): Seq[DBIO[Unit]] = {
    keys.toSeq map toScopeAction(workflowExecution, backend)
  }

  override def insertCalls(workflowId: WorkflowId, keys: Traversable[ExecutionStoreKey], backend: Backend)
                          (implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      workflowExecution <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.head
      _ <- DBIO.sequence(toScopeActions(workflowExecution, backend, keys))
    } yield ()

    runTransaction(action)
  }

  private def insertEmptyExecutionInfos(insertResult: this.dataAccess.executionsAutoInc.SingleInsertResult,
                                   backend: Backend): DBIO[_] = {

    def buildInsert(key: String) = {
      dataAccess.executionInfosAutoInc += new ExecutionInfo(insertResult.executionId.get, key, None)
    }

    val inserts = backend.executionInfoKeys map buildInsert
    DBIO.sequence(inserts)
  }

  // Converts a single Call to a composite DBIOAction[] that inserts the correct rows
  private def toScopeAction(workflowExecution: WorkflowExecution, backend: Backend)
                           (key: ExecutionStoreKey)(implicit ec: ExecutionContext): DBIO[Unit] = {
    for {
      // Insert an execution row
      executionInsert <- dataAccess.executionsAutoInc +=
        new Execution(
          workflowExecutionId = workflowExecution.workflowExecutionId.get,
          callFqn = key.scope.fullyQualifiedName,
          index = key.index.fromIndex,
          status = ExecutionStatus.NotStarted.toString,
          rc = None,
          startDt = None,
          endDt = None,
          backendType = backend.backendType.displayName,
          attempt = key.attempt)

      // Add the empty execution info rows using the execution as the FK.
      _ <- insertEmptyExecutionInfos(executionInsert, backend)
    } yield ()
  }

  override def getWorkflowState(workflowId: WorkflowId)
                               (implicit ec: ExecutionContext): Future[Option[WorkflowState]] = {
    val action = for {
      maybeWorkflowExecution <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(
        workflowId.id.toString).result.headOption
      workflowState = maybeWorkflowExecution map { w => WorkflowState.fromString(w.status) }
    } yield workflowState

    runTransaction(action)
  }

  override def getExecutionStatuses(workflowId: WorkflowId)
                                   (implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, CallStatus]] = {
    val action = for {
    // NOTE: For now, intentionally causes query to error out instead of returning an Map.empty
      workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(
        workflowId.id.toString).result.head

      // Alternatively, could use a dataAccess.executionCallFqnsAndStatusesByWorkflowExecutionUuid
      executionKeyAndStatusResults <- dataAccess.executionCallFqnsAndStatusesByWorkflowExecutionId(
        workflowExecutionResult.workflowExecutionId.get).result

      executionStatuses = executionKeyAndStatusResults map { execution =>
        (ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex, execution.attempt),
         CallStatus(execution.status.toExecutionStatus, execution.rc, execution.overallHash map { ExecutionHash(_, execution.dockerImageHash )}, None)) }

    } yield executionStatuses.toMap

    runTransaction(action)
  }

  override def getExecutionStatuses(workflowId: WorkflowId, fqn: FullyQualifiedName)
                                   (implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, CallStatus]] = {
    val action = for {
      workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(
        workflowId.toString).result.head

      executionKeyAndStatusResults <- dataAccess.executionStatusByWorkflowExecutionIdAndCallFqn(
        (workflowExecutionResult.workflowExecutionId.get, fqn)).result

      executionStatuses = executionKeyAndStatusResults map { execution =>
        (ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex, execution.attempt), CallStatus(execution.status.toExecutionStatus, execution.rc,
          execution.overallHash map { ExecutionHash(_, execution.dockerImageHash) }, None)) }
    } yield executionStatuses.toMap

    runTransaction(action)
  }

  override def getExecutionStatus(workflowId: WorkflowId, key: ExecutionDatabaseKey)
                                 (implicit ec: ExecutionContext): Future[Option[CallStatus]] = {
    val action = for {
      workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.head

      executionStatuses <- dataAccess.executionStatusesAndReturnCodesByWorkflowExecutionIdAndCallKey(
        (workflowExecutionResult.workflowExecutionId.get, key.fqn, key.index.fromIndex, key.attempt)).result

      maybeStatus = executionStatuses.headOption map { case (status, rc, hash, dockerHash) => CallStatus(status.toExecutionStatus, rc, hash map { ExecutionHash(_, dockerHash )}, None) }
    } yield maybeStatus
    runTransaction(action)
  }

  override def getWorkflow(workflowExecutionId: Int)(implicit ec: ExecutionContext): Future[WorkflowDescriptor] = {
    val action = for {
      workflowExecutionResult <- dataAccess.workflowExecutionsByPrimaryKey(workflowExecutionId).result.head
      workflowAux <- dataAccess.workflowExecutionAuxesByWorkflowExecutionUuid(workflowExecutionResult.workflowExecutionUuid).result.head
      workflowDescriptor = WorkflowDescriptor(
        WorkflowId(UUID.fromString(workflowExecutionResult.workflowExecutionUuid)),
        WorkflowSourceFiles(workflowAux.wdlSource.toRawString, workflowAux.jsonInputs.toRawString, workflowAux.workflowOptions.toRawString)
      )
    } yield workflowDescriptor

    runTransaction(action)
  }

  override def getWorkflow(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowDescriptor] = {
    val action = for {
      workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.head
      workflowAux <- dataAccess.workflowExecutionAuxesByWorkflowExecutionUuid(workflowExecutionResult.workflowExecutionUuid).result.head
      workflowDescriptor = WorkflowDescriptor(
        workflowId,
        WorkflowSourceFiles(workflowAux.wdlSource.toRawString, workflowAux.jsonInputs.toRawString, workflowAux.workflowOptions.toRawString)
      )
    } yield workflowDescriptor

    runTransaction(action)
  }

  override def getWorkflowsByState(states: Traversable[WorkflowState])
                                  (implicit ec: ExecutionContext): Future[Traversable[WorkflowDescriptor]] = {
    val action = for {
      workflowExecutionResults <- dataAccess.workflowExecutionsByStatuses(states.map(_.toString)).result

      workflowDescriptors <- DBIO.sequence(
        workflowExecutionResults map { workflowExecutionResult =>
          val workflowExecutionAuxResult = dataAccess.workflowExecutionAuxesByWorkflowExecutionUuid(
            workflowExecutionResult.workflowExecutionUuid).result.head

          workflowExecutionAuxResult map { workflowExecutionAux =>
            WorkflowDescriptor(
              WorkflowId.fromString(workflowExecutionResult.workflowExecutionUuid),
              WorkflowSourceFiles(
                workflowExecutionAux.wdlSource.toRawString,
                workflowExecutionAux.jsonInputs.toRawString,
                workflowExecutionAux.workflowOptions.toRawString
              )
            )
          }
        }
      )

    } yield workflowDescriptors

    runTransaction(action)
  }

  override def updateWorkflowState(workflowId: WorkflowId, workflowState: WorkflowState)
                                  (implicit ec: ExecutionContext): Future[Unit] = {
    val query = dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).extract
    val endDate = if (workflowState.isTerminal) Option(new Date().toTimestamp) else None

    val action = for {
      count <- query.map(w => (w.status, w.endDt)).update(workflowState.toString, endDate)
      _ = require(count == 1, s"Unexpected workflow execution update count $count")
    } yield ()

    runTransaction(action)
  }

  override def getExecutionInfos(workflowId: WorkflowId, call: Call, attempt: Int)(implicit ec: ExecutionContext): Future[Traversable[ExecutionInfo]] = {
    val action = for {
      executionResult <- dataAccess.executionsByWorkflowExecutionUuidAndCallFqnAndAttempt(
        (workflowId.toString, call.fullyQualifiedName, attempt)).result.head

      executionInfos <- dataAccess.executionInfosByExecutionId(executionResult.executionId.get).result

    } yield executionInfos

    runTransaction(action)
  }

  override def updateExecutionInfo(workflowId: WorkflowId,
                                   callKey: BackendCallKey,
                                   key: String,
                                   value: Option[String])(implicit ec: ExecutionContext): Future[Unit] = {
    import ExecutionIndex._
    val action = for {
      executionResult <- dataAccess.executionsByWorkflowExecutionUuidAndCallFqnAndShardIndexAndAttempt(
        workflowId.toString, callKey.scope.fullyQualifiedName, callKey.index.fromIndex, callKey.attempt).result.head

      backendUpdate <- dataAccess.executionInfoValueByExecutionAndKey(executionResult.executionId.get, key).update(value)
      _ = require(backendUpdate == 1, s"Unexpected backend update count $backendUpdate")
    } yield ()

    runTransaction(action)
  }

  private def toSymbolStoreEntries(symbolResults: Traversable[Symbol]) =
    symbolResults map toSymbolStoreEntry

  private def toSymbolStoreEntry(symbolResult: Symbol) = {
    val wdlType = WdlType.fromWdlString(symbolResult.wdlType)
    new SymbolStoreEntry(
      new SymbolStoreKey(
        symbolResult.scope,
        symbolResult.name,
        symbolResult.index.toIndex,
        input = symbolResult.io == IoInput // input = true, if db contains "INPUT"
      ),
      wdlType,
      symbolResult.wdlValue map { v => dbEntryToWdlValue(v.toRawString, wdlType) },
      symbolResult.symbolHash map SymbolHash
    )
  }

  override def getAllSymbolStoreEntries(workflowId: WorkflowId)
                                       (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    val action = dataAccess.allSymbols(workflowId.toString).result
    runTransaction(action) map toSymbolStoreEntries
  }

  /** Get all inputs for the scope of this key. */
  override def getInputs(workflowId: WorkflowId, call: Call)
                        (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    require(call != null, "call cannot be null")
    getSymbols(workflowId, IoInput, Option(call.fullyQualifiedName))
  }

  /** Get all outputs for the scope of this key. */
  override def getOutputs(workflowId: WorkflowId, key: ExecutionDatabaseKey)
                         (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    require(key != null, "key cannot be null")
    getSymbols(workflowId, IoOutput, Option(key.fqn), key.index)
  }

  /** Returns all NON SHARDS outputs for this workflowId */
  override def getWorkflowOutputs(workflowId: WorkflowId)
                                 (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    val action = dataAccess.symbolsForWorkflowOutput(workflowId.toString).result
    runTransaction(action) map toSymbolStoreEntries
  }

  private def getSymbols(workflowId: WorkflowId,
                         ioValue: IoValue,
                         callFqnOption: Option[FullyQualifiedName] = None,
                         callIndexOption: Option[Int] = None)
                        (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    val action = dataAccess.symbolsByWorkflowExecutionUuidAndIoAndMaybeScope(
      workflowId.toString, ioValue, callFqnOption, callIndexOption
    ).result

    runTransaction(action) map toSymbolStoreEntries
  }

  /** Should fail if a value is already set.  The keys in the Map are locally qualified names. */
  override def setOutputs(workflowId: WorkflowId, key: OutputKey, callOutputs: WorkflowOutputs,
                          reportableResults: Seq[ReportableSymbol])
                         (implicit ec: ExecutionContext): Future[Unit] = {
    val reportableResultNames = reportableResults map { _.fullyQualifiedName }
    val action = for {
      workflowExecution <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.head
      _ <- dataAccess.symbolsAutoInc ++= callOutputs map {
        case (symbolLocallyQualifiedName, CallOutput(wdlValue, hash)) =>
          val reportableSymbol = key.index.fromIndex == -1 && reportableResultNames.contains(key.scope.fullyQualifiedName + "." + symbolLocallyQualifiedName)
          new Symbol(
            workflowExecution.workflowExecutionId.get,
            key.scope.fullyQualifiedName,
            symbolLocallyQualifiedName,
            key.index.fromIndex,
            IoOutput,
            reportableSymbol,
            wdlValue.wdlType.toWdlString,
            Option(wdlValueToDbValue(wdlValue).toClob),
            hash.value.map(_.value)
          )
      }
    } yield ()

    runTransaction(action)
  }

  /**
    * Updates the existing input symbols to replace expressions with real values.
    *
    * @return The number of rows updated - as a Future.
    */
  override def updateCallInputs(workflowId: WorkflowId, key: BackendCallKey, callInputs: CallInputs)
                               (implicit ec: ExecutionContext): Future[Int] = {
    type ProjectionFunction = SlickDataAccess.this.dataAccess.Symbols => (Rep[String], Rep[Option[Clob]])
    val projectionFn: ProjectionFunction = (s: SlickDataAccess.this.dataAccess.Symbols) => (s.wdlType, s.wdlValue)

    val inputUpdateActions = callInputs map { case (inputName, wdlValue) =>
      for {
        workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.head
        symbols = dataAccess.symbolsFilterByWorkflowAndScopeAndNameAndIndex(workflowExecutionResult.workflowExecutionId.get, key.scope.fullyQualifiedName, inputName, key.index.fromIndex)
        count <- symbols.map(projectionFn).update(wdlValue.wdlType.toWdlString, Option(wdlValueToDbValue(wdlValue).toClob))
      } yield count
    }

    // Do an FP dance to get the DBIOAction[Iterable[Int]] from Iterable[DBIOAction[Int]].
    val allInputUpdatesAction = DBIO.sequence(inputUpdateActions)
    runTransaction(allInputUpdatesAction) map { _.sum }
  }

  override def setExecutionEvents(workflowId: WorkflowId, callFqn: String, shardIndex: Option[Int], attempt: Int,
                                  events: Seq[ExecutionEventEntry])(implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      execution <- shardIndex match {
        case Some(idx) => dataAccess.executionsByWorkflowExecutionUuidAndCallFqnAndShardIndexAndAttempt(workflowId.toString, callFqn, idx, attempt).result.head
        case None => dataAccess.executionsByWorkflowExecutionUuidAndCallFqnAndAttempt(workflowId.toString, callFqn, attempt).result.head
      }
      _ <- dataAccess.executionEventsAutoInc ++= events map { executionEventEntry =>
        new ExecutionEvent(
          execution.executionId.get,
          executionEventEntry.description,
          new Timestamp(executionEventEntry.startTime.getMillis),
          new Timestamp(executionEventEntry.endTime.getMillis))
      }
    } yield ()

    runTransaction(action)
  }

  override def getAllExecutionEvents(workflowId: WorkflowId)(implicit ec: ExecutionContext):
  Future[Map[ExecutionDatabaseKey, Seq[ExecutionEventEntry]]] = {
    // The database query gives us a Seq[(CallFqn, ExecutionEvent)]. We want a Map[CallFqn -> ExecutionEventEntry].
    // So let's do some functional programming!
    val action = dataAccess.executionEventsByWorkflowExecutionUuid(workflowId.toString).result
    runTransaction(action) map toExecutionEvents
  }

  private def toExecutionEvents(events: Traversable[((String, Int, Int), ExecutionEvent)])
                               (implicit ec: ExecutionContext): Map[ExecutionDatabaseKey, Seq[ExecutionEventEntry]] = {
      // First: Group all the entries together by name
      val grouped: Map[ExecutionDatabaseKey, Seq[((String, Int, Int), ExecutionEvent)]] = events.toSeq groupBy { case ((fqn: String, idx: Int, attempt: Int), event: ExecutionEvent) => ExecutionDatabaseKey(fqn, idx.toIndex, attempt) }
      // Second: Transform the values. The value no longer needs the String since that's now part of the Map, and
      // convert the executionEvent into a friendlier ExecutionEventEntry:
      grouped mapValues { _ map { case (_ , event: ExecutionEvent) =>
        ExecutionEventEntry(
          event.description,
          new DateTime(event.startTime.getTime),
          new DateTime(event.endTime.getTime))
      } }
  }

  private def executionsForStatusUpdate(workflowId: WorkflowId, scopeKeys: Traversable[ExecutionDatabaseKey])(implicit ec: ExecutionContext) = for {
    workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.head
    executions = dataAccess.executionsByWorkflowExecutionIdAndScopeKeys(workflowExecutionResult.workflowExecutionId.get, scopeKeys)
  } yield executions

  private def assertUpdateCount(updates: Int, expected: Int) = require(updates == expected, s"Execution update count $updates did not match scopes size $expected")

  def setStartingStatus(workflowId: WorkflowId, scopeKeys: Traversable[ExecutionDatabaseKey])(implicit ec: ExecutionContext): Future[Unit] = {
    if (scopeKeys.isEmpty) Future.successful(()) else runTransaction(setStartingStatusAction(workflowId, scopeKeys))
  }

  def setStartingStatusAction(workflowId: WorkflowId, scopeKeys: Traversable[ExecutionDatabaseKey])(implicit ec: ExecutionContext) = {
    for {
      executions <- executionsForStatusUpdate(workflowId, scopeKeys)
      count <- executions.map(e => (e.status, e.startDt)).update((ExecutionStatus.Starting.toString, Option(new Date().toTimestamp)))
      _ = assertUpdateCount(count, scopeKeys.size)
    } yield ()
  }

  def updateStatus(workflowId: WorkflowId, scopeKeys: Traversable[ExecutionDatabaseKey], status: ExecutionStatus)(implicit ec: ExecutionContext): Future[Unit] = {
    if (scopeKeys.isEmpty) Future.successful(()) else runTransaction(updateStatusAction(workflowId, scopeKeys, status))
  }

  def updateStatusAction(workflowId: WorkflowId, scopeKeys: Traversable[ExecutionDatabaseKey], status: ExecutionStatus)(implicit ec: ExecutionContext) = {
    require(!status.isTerminal && !(status == Starting))

    for {
      executions <- executionsForStatusUpdate(workflowId, scopeKeys)
      count <- executions.map(_.status).update(status.toString)
      _ = assertUpdateCount(count, scopeKeys.size)
    } yield ()
  }

  def setTerminalStatus(workflowId: WorkflowId, scopeKey: ExecutionDatabaseKey, status: ExecutionStatus,
                        scriptReturnCode: Option[Int], hash: Option[ExecutionHash],
                        resultsClonedFrom: Option[BackendCall])(implicit ec: ExecutionContext): Future[Unit] = {
    runTransaction(setTerminalStatusAction(workflowId, scopeKey, status, scriptReturnCode, hash, resultsClonedFrom))
  }

  def setTerminalStatusAction(workflowId: WorkflowId, scopeKey: ExecutionDatabaseKey, status: ExecutionStatus,
                              scriptReturnCode: Option[Int], hash: Option[ExecutionHash],
                              resultsClonedFrom: Option[BackendCall])(implicit ec: ExecutionContext) = {
    require(status.isTerminal)

    val overallHash = hash map { _.overallHash }
    val dockerHash = hash flatMap { _.dockerHash }
    val projection = { e: SlickDataAccess.this.dataAccess.Executions => (e.status, e.endDt, e.rc, e.executionHash, e.dockerImageHash, e.resultsClonedFrom) }

    val findResultsClonedFromId = resultsClonedFrom map { backendCall =>
      for {
        workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(backendCall.workflowDescriptor.id.toString).result.head
        execution <- dataAccess.executionsByWorkflowExecutionIdAndCallFqnAndIndexAndAttempt(
          workflowExecutionResult.workflowExecutionId.get, backendCall.key.scope.fullyQualifiedName, backendCall.key.index.fromIndex, backendCall.key.attempt).result.head
      } yield Option(execution.executionId.get)
    } getOrElse DBIO.successful(None)

    for {
      executions <- executionsForStatusUpdate(workflowId, List(scopeKey))
      clonedFromId <- findResultsClonedFromId
      count <- executions.map(projection).update((status.toString, Option(new Date().toTimestamp), scriptReturnCode, overallHash, dockerHash, clonedFromId))
      _ = assertUpdateCount(count, 1)
    } yield ()
  }

  override def getExecutions(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[Execution]] = {
    val action = dataAccess.executionsByWorkflowExecutionUuid(id.toString).result

    runTransaction(action)
  }

  override def getExecutionsForRestart(id: WorkflowId)
                                      (implicit ec: ExecutionContext): Future[Traversable[Execution]] = {
    val action = dataAccess.executionsForRestartByWorkflowExecutionUuid(id.toString).result

    runTransaction(action)
  }

  override def getExecutionsWithResuableResultsByHash(hash: String)
                                                     (implicit ec: ExecutionContext): Future[Traversable[Execution]] = {
    val action = dataAccess.executionsWithReusableResultsByExecutionHash(hash).result

    runTransaction(action)
  }

  override def getWorkflowExecution(workflowId: WorkflowId)
                                   (implicit ec: ExecutionContext): Future[WorkflowExecution] = {
    val action = dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.headOption

    runTransaction(action) map { _.getOrElse(throw new NoSuchElementException(s"Workflow $workflowId not found.")) }
  }

  override def getWorkflowExecutionAux(id: WorkflowId)
                                      (implicit ec: ExecutionContext): Future[WorkflowExecutionAux] = {
    val action = dataAccess.workflowExecutionAuxesByWorkflowExecutionUuid(id.toString).result.headOption

    runTransaction(action) map { _.getOrElse(throw new NoSuchElementException(s"No workflow execution aux found for ID '$id'.")) }
  }

  override def getAllInputs(workflowId: WorkflowId)
                           (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    val action = dataAccess.symbolsByWorkflowExecutionUuidAndIo(workflowId.toString, IoInput).result

    runTransaction(action) map toSymbolStoreEntries
  }

  override def getAllOutputs(workflowId: WorkflowId)
                            (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]] = {
    val action = dataAccess.symbolsByWorkflowExecutionUuidAndIo(workflowId.toString, IoOutput).result

    runTransaction(action) map toSymbolStoreEntries
  }

  override def updateWorkflowOptions(workflowId: WorkflowId, workflowOptionsJson: String)(implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      workflowExecution <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.id.toString).result.head
      count <- dataAccess.workflowOptionsFromWorkflowId(workflowExecution.workflowExecutionId.get).update(workflowOptionsJson.toClob)
      _ = require(count == 1, s"Unexpected workflow aux update count $count")
    } yield ()

    runTransaction(action)
  }

  override def resetTransientExecutions(workflowId: WorkflowId, isResetable: (Execution, Seq[ExecutionInfo]) => Boolean)(implicit ec: ExecutionContext): Future[Unit] = {
    val action = for {
      tuples <- dataAccess.executionsAndExecutionInfosByWorkflowId(workflowId.toString).result
      executionsAndInfos = tuples groupBy { _._1 } mapValues { _ map { _._2 } }

      transientDatabaseKeys = (executionsAndInfos filter Function.tupled(isResetable)).keys map { _.toKey }
      _ <- updateStatusAction(workflowId, transientDatabaseKeys, ExecutionStatus.NotStarted)
    } yield ()

    runTransaction(action)
  }

  override def findResumableExecutions(workflowId: WorkflowId,
                                       isResumable: (Execution, Seq[ExecutionInfo]) => Boolean,
                                       jobKeyBuilder: (Execution, Seq[ExecutionInfo]) => JobKey)
                                      (implicit ec: ExecutionContext): Future[Traversable[ExecutionKeyToJobKey]] = {
    val action = for {
      tuples <- dataAccess.executionsAndExecutionInfosByWorkflowId(workflowId.toString).result
      executionsAndInfos = tuples groupBy { _._1 } mapValues { _ map { _._2 } }

      resumablePairs = executionsAndInfos collect { case (e, ei) if isResumable(e, ei) => e.toKey -> jobKeyBuilder(e, ei) }
    } yield resumablePairs

    runTransaction(action) map { _ map { case (e, j) => ExecutionKeyToJobKey(e, j) } toTraversable }
  }

  override def queryWorkflows(queryParameters: WorkflowQueryParameters)
                             (implicit ec: ExecutionContext): Future[WorkflowQueryResponse] = {
    val action = dataAccess.queryWorkflowExecutions(queryParameters).result
    runTransaction(action) map { workflows =>
      WorkflowQueryResponse(workflows map { workflow =>
        WorkflowQueryResult(
          id = workflow.workflowExecutionUuid,
          name = workflow.name,
          status = workflow.status,
          start = new DateTime(workflow.startDt),
          end = workflow.endDt map { new DateTime(_) })
      })
    }
  }

  override def updateCallCaching(parameters: CallCachingParameters)(implicit ec: ExecutionContext): Future[Int] = {
    // Figure out which of the three possible queries to use based on whether a call has been specified and
    // if so whether an index has been specified.
    val executionQuery: (Int) => Query[dataAccess.Executions, Execution, Seq] = {
      (parameters.callKey, parameters.callKey flatMap { _.index }) match {
        case (Some(key), Some(idx)) => dataAccess.executionsByWorkflowExecutionIdAndCallFqnAndIndexAndAttempt(_: Int, key.fqn, idx, key.attempt).extract
        case (Some(key), None) => dataAccess.executionsByWorkflowExecutionIdAndCallFqnAndAttempt(_: Int, key.fqn, key.attempt).extract
        case _ => dataAccess.executionsByWorkflowExecutionId(_: Int).extract
      }
    }

    val action = for {
      workflowExecution <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(parameters.workflowId.id.toString).result.head
      count <- executionQuery(workflowExecution.workflowExecutionId.get).map(_.allowsResultReuse).update(parameters.allow)
    } yield count

    runTransaction(action)
  }

  override def infosByExecution(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[ExecutionInfosByExecution]] = {
    val action = for {
      rawTuples <- dataAccess.executionsAndExecutionInfosByWorkflowId(id.toString).result
      // Group by execution, then remove the execution keys of the execution -> execution info tuple.
      // The net result is Execution -> Seq[ExecutionInfo].
      groupedTuples = rawTuples groupBy { _._1 } mapValues { _ map { _._2 } }
      infosByExecution = groupedTuples map ExecutionInfosByExecution.tupled
    } yield infosByExecution

    runTransaction(action)
  }
}
