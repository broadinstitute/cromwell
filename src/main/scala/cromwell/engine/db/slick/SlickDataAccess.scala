package cromwell.engine.db.slick

import java.sql.{Clob, Timestamp}
import java.util.{Date, UUID}
import javax.sql.rowset.serial.SerialClob

import _root_.slick.backend.DatabaseConfig
import _root_.slick.driver.JdbcProfile
import com.typesafe.config.{Config, ConfigFactory, ConfigValueFactory}
import cromwell.engine.ExecutionIndex._
import cromwell.engine.ExecutionStatus._
import cromwell.engine.backend.jes.{JesBackend, JesJobKey}
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.backend.sge.SgeBackend
import cromwell.engine.backend.{Backend, WorkflowQueryResult}
import cromwell.engine.db._
import cromwell.engine.finalcall.FinalCall
import cromwell.engine.workflow._
import cromwell.engine.{WorkflowOutputs, _}
import cromwell.webservice.{CallCachingParameters, WorkflowQueryParameters, WorkflowQueryResponse}
import lenthall.config.ScalaConfig._
import org.joda.time.DateTime
import org.slf4j.LoggerFactory
import wdl4s.types.{WdlPrimitiveType, WdlType}
import wdl4s.values.WdlValue
import wdl4s.{CallInputs, _}

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

  // NOTE: Used for slick flatMap. May switch to custom ExecutionContext the future
  private implicit val executionContext = ExecutionContext.global

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

  override def close(): Unit = database.close()

  // Run action with an outer transaction
  private def runTransaction[R](action: DBIOAction[R, _ <: NoStream, _ <: Effect]): Future[R] = {
    database.run(action.transactionally)
  }

  /**
   * Creates a row in each of the backend-info specific tables for each key in `keys` corresponding to the backend
   * `backend`.  Or perhaps defer this?
   */
  override def createWorkflow(workflowDescriptor: WorkflowDescriptor,
                              workflowInputs: Traversable[SymbolStoreEntry],
                              scopes: Traversable[Scope],
                              backend: Backend): Future[Unit] = {

    val scopeKeys: Traversable[ExecutionStoreKey] = scopes collect {
      case call: Call => BackendCallKey(call, None)
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
                            keys: Traversable[ExecutionStoreKey]): Seq[DBIO[Unit]] = {
    keys.toSeq map toScopeAction(workflowExecution, backend)
  }

  override def insertCalls(workflowId: WorkflowId, keys: Traversable[ExecutionStoreKey], backend: Backend): Future[Unit] = {
    val action = for {
      workflowExecution <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.head
      _ <- DBIO.sequence(toScopeActions(workflowExecution, backend, keys))
    } yield ()

    runTransaction(action)
  }

  // Converts a single Call to a composite DBIOAction[] that inserts the correct rows
  private def toScopeAction(workflowExecution: WorkflowExecution, backend: Backend)
                           (key: ExecutionStoreKey): DBIO[Unit] = {
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
          endDt = None)

      // Depending on the backend, insert a job specific row
      _ <- backend match {
        case _: LocalBackend =>
          dataAccess.localJobsAutoInc +=
            new LocalJob(
              executionInsert.executionId.get,
              None,
              None)
        case j: JesBackend =>
          // FIXME: Placeholder for now, discussed w/ Khalid
          dataAccess.jesJobsAutoInc += new JesJob(executionInsert.executionId.get, None, None, None)
        case s: SgeBackend =>
          dataAccess.sgeJobsAutoInc += new SgeJob(executionInsert.executionId.get, None)
        case null =>
          throw new IllegalArgumentException("Backend is null")
        case unknown =>
          throw new IllegalArgumentException("Unknown backend: " + backend.getClass)
      }

    } yield ()
  }

  override def getWorkflowState(workflowId: WorkflowId): Future[Option[WorkflowState]] = {
    val action = for {
      maybeWorkflowExecution <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(
        workflowId.id.toString).result.headOption
      workflowState = maybeWorkflowExecution map { w => WorkflowState.fromString(w.status) }
    } yield workflowState

    runTransaction(action)
  }

  override def getExecutionStatuses(workflowId: WorkflowId): Future[Map[ExecutionDatabaseKey, CallStatus]] = {
    val action = for {
    // NOTE: For now, intentionally causes query to error out instead of returning an Map.empty
      workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(
        workflowId.id.toString).result.head

      // Alternatively, could use a dataAccess.executionCallFqnsAndStatusesByWorkflowExecutionUuid
      executionKeyAndStatusResults <- dataAccess.executionCallFqnsAndStatusesByWorkflowExecutionId(
        workflowExecutionResult.workflowExecutionId.get).result

      executionStatuses = executionKeyAndStatusResults map {
        case (fqn, indexInt, status, rc, executionHash, dockerHash) => (ExecutionDatabaseKey(fqn, indexInt.toIndex), CallStatus(status.toExecutionStatus, rc, executionHash map { ExecutionHash(_, dockerHash )}, None)) }

    } yield executionStatuses.toMap

    runTransaction(action)
  }

  override def getExecutionStatuses(workflowId: WorkflowId, fqn: FullyQualifiedName): Future[Map[ExecutionDatabaseKey, CallStatus]] = {
    val action = for {
      workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(
        workflowId.toString).result.head

      executionKeyAndStatusResults <- dataAccess.executionStatusByWorkflowExecutionIdAndCallFqn(
        (workflowExecutionResult.workflowExecutionId.get, fqn)).result

      executionStatuses = executionKeyAndStatusResults map { case (callFqn, indexInt, status, rc, hash, dockerHash) =>
        (ExecutionDatabaseKey(callFqn, indexInt.toIndex), CallStatus(status.toExecutionStatus, rc,
          hash map { ExecutionHash(_, dockerHash) }, None)) }
    } yield executionStatuses.toMap

    runTransaction(action)
  }

  override def getExecutionStatus(workflowId: WorkflowId, key: ExecutionDatabaseKey): Future[Option[CallStatus]] = {
    val action = for {
      workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.head

      executionStatuses <- dataAccess.executionStatusesAndReturnCodesByWorkflowExecutionIdAndCallKey(
        (workflowExecutionResult.workflowExecutionId.get, key.fqn, key.index.fromIndex)).result

      maybeStatus = executionStatuses.headOption map { case (status, rc, hash, dockerHash) => CallStatus(status.toExecutionStatus, rc, hash map { ExecutionHash(_, dockerHash )}, None) }
    } yield maybeStatus
    runTransaction(action)
  }

  override def getWorkflow(workflowExecutionId: Int): Future[WorkflowDescriptor] = {
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

  override def getWorkflow(workflowId: WorkflowId): Future[WorkflowDescriptor] = {
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

  override def getWorkflowsByState(states: Traversable[WorkflowState]): Future[Traversable[WorkflowDescriptor]] = {
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

  override def updateWorkflowState(workflowId: WorkflowId, workflowState: WorkflowState): Future[Unit] = {
    val query = dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).extract
    val endDate = if (workflowState.isTerminal) Option(new Date().toTimestamp) else None

    val action = for {
      count <- query.map(w => (w.status, w.endDt)).update(workflowState.toString, endDate)
      _ = require(count == 1, s"Unexpected workflow execution update count $count")
    } yield ()

    runTransaction(action)
  }

  override def getExecutionBackendInfo(workflowId: WorkflowId, call: Call): Future[CallBackendInfo] = {
    val action = for {
      executionResult <- dataAccess.executionsByWorkflowExecutionUuidAndCallFqn(
        (workflowId.toString, call.fullyQualifiedName)).result.head

      localJobResultOption <- dataAccess.localJobsByExecutionId(executionResult.executionId.get).result.headOption

      jesJobResultOption <- dataAccess.jesJobsByExecutionId(executionResult.executionId.get).result.headOption

      sgeJobResultOption <- dataAccess.sgeJobsByExecutionId(executionResult.executionId.get).result.headOption

      jobResultOption = localJobResultOption orElse jesJobResultOption orElse sgeJobResultOption
      backendInfo = jobResultOption match {
        case Some(localJobResult: LocalJob) =>
          new LocalCallBackendInfo(localJobResult.pid)
        case Some(jesJobResult: JesJob) =>
          new JesCallBackendInfo(jesJobResult.jesRunId map JesId,
            jesJobResult.jesStatus map JesStatus)
        case Some(sgeJobResult: SgeJob) =>
          new SgeCallBackendInfo(sgeJobResult.sgeJobNumber)
        case _ =>
          throw new IllegalArgumentException(
            s"Unknown backend from db for (uuid, fqn): " +
              s"($workflowId, ${call.fullyQualifiedName})")
      }

    } yield backendInfo

    runTransaction(action)
  }

  override def updateExecutionBackendInfo(workflowId: WorkflowId,
                                          callKey: BackendCallKey,
                                          backendInfo: CallBackendInfo): Future[Unit] = {
    require(backendInfo != null, "backend info is null")

    import ExecutionIndex._
    val action = for {
      executionResult <- dataAccess.executionsByWorkflowExecutionUuidAndCallFqnAndShardIndex(
        workflowId.toString, callKey.scope.fullyQualifiedName, callKey.index.fromIndex).result.head

      backendUpdate <- backendInfo match {
        case localBackendInfo: LocalCallBackendInfo =>
          dataAccess.localJobPidsByExecutionId(
            executionResult.executionId.get).update(
              localBackendInfo.processId)

        case jesBackendInfo: JesCallBackendInfo =>
          dataAccess.jesIdsAndJesStatusesByExecutionId(
            executionResult.executionId.get).update(
              jesBackendInfo.jesId map { _.id },
              jesBackendInfo.jesStatus map { _.status }
            )

        case sgeBackendInfo: SgeCallBackendInfo =>
          dataAccess.sgeJobNumberByExecutionId(
            executionResult.executionId.get
          ).update(sgeBackendInfo.sgeJobNumber)
      }

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

  override def getAllSymbolStoreEntries(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]] = {
    val action = dataAccess.allSymbols(workflowId.toString).result
    runTransaction(action) map toSymbolStoreEntries
  }

  /** Get all inputs for the scope of this key. */
  override def getInputs(workflowId: WorkflowId, call: Call): Future[Traversable[SymbolStoreEntry]] = {
    require(call != null, "call cannot be null")
    getSymbols(workflowId, IoInput, Option(call.fullyQualifiedName))
  }

  /** Get all outputs for the scope of this key. */
  override def getOutputs(workflowId: WorkflowId, key: ExecutionDatabaseKey): Future[Traversable[SymbolStoreEntry]] = {
    require(key != null, "key cannot be null")
    getSymbols(workflowId, IoOutput, Option(key.fqn), key.index)
  }

  /** Returns all NON SHARDS outputs for this workflowId */
  override def getWorkflowOutputs(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]] = {
    val action = dataAccess.symbolsForWorkflowOutput(workflowId.toString).result
    runTransaction(action) map toSymbolStoreEntries
  }

  private def getSymbols(workflowId: WorkflowId,
                         ioValue: IoValue,
                         callFqnOption: Option[FullyQualifiedName] = None,
                         callIndexOption: Option[Int] = None): Future[Traversable[SymbolStoreEntry]] = {
    val action = dataAccess.symbolsByWorkflowExecutionUuidAndIoAndMaybeScope(
      workflowId.toString, ioValue, callFqnOption, callIndexOption
    ).result

    runTransaction(action) map toSymbolStoreEntries
  }

  /** Should fail if a value is already set.  The keys in the Map are locally qualified names. */
  override def setOutputs(workflowId: WorkflowId, key: OutputKey, callOutputs: WorkflowOutputs, reportableResults: Seq[ReportableSymbol]): Future[Unit] = {
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
    * @return The number of rows updated - as a Future.
    */
  override def updateCallInputs(workflowId: WorkflowId, key: BackendCallKey, callInputs: CallInputs): Future[Int] = {
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

  override def setExecutionEvents(workflowId: WorkflowId, callFqn: String, shardIndex: Option[Int], events: Seq[ExecutionEventEntry]): Future[Unit] = {
    val action = for {
      execution <- shardIndex match {
        case Some(idx) => dataAccess.executionsByWorkflowExecutionUuidAndCallFqnAndShardIndex(workflowId.toString, callFqn, idx).result.head
        case None => dataAccess.executionsByWorkflowExecutionUuidAndCallFqn(workflowId.toString, callFqn).result.head
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

  override def getAllExecutionEvents(workflowId: WorkflowId): Future[Map[ExecutionDatabaseKey, Seq[ExecutionEventEntry]]] = {
    // The database query gives us a Seq[(CallFqn, ExecutionEvent)]. We want a Map[CallFqn -> ExecutionEventEntry].
    // So let's do some functional programming!
    val action = dataAccess.executionEventsByWorkflowExecutionUuid(workflowId.toString).result
    runTransaction(action) map toExecutionEvents
  }

  private def toExecutionEvents(events: Traversable[((String, Int), ExecutionEvent)]): Map[ExecutionDatabaseKey, Seq[ExecutionEventEntry]] = {
      // First: Group all the entries together by name
      val grouped: Map[ExecutionDatabaseKey, Seq[((String, Int), ExecutionEvent)]] = events.toSeq groupBy { case ((fqn: String, idx: Int), event: ExecutionEvent) => ExecutionDatabaseKey(fqn, idx.toIndex) }
      // Second: Transform the values. The value no longer needs the String since that's now part of the Map, and
      // convert the executionEvent into a friendlier ExecutionEventEntry:
      grouped mapValues { _ map { case (_ , event: ExecutionEvent) =>
        ExecutionEventEntry(
          event.description,
          new DateTime(event.startTime.getTime),
          new DateTime(event.endTime.getTime))
      } }
  }

  private def setStatusAction(workflowId: WorkflowId, scopeKeys: Traversable[ExecutionDatabaseKey],
                              callStatus: CallStatus): DBIO[Unit] = {
    // Describes a function from an input `Executions` to a projection of fields to be updated.
    type ProjectionFunction = SlickDataAccess.this.dataAccess.Executions => (
      Rep[String], Rep[Option[Timestamp]], Rep[Option[Int]], Rep[Option[String]], Rep[Option[String]], Rep[Option[Int]])
    // If the call status is Starting, target the start date for update, otherwise target the end date.  The end date
    // is only set to a non-None value if the status is terminal.
    val projectionFn: ProjectionFunction = if (callStatus.isStarting)
      e => (e.status, e.startDt, e.rc, e.executionHash, e.dockerImageHash, e.resultsClonedFrom)
    else
      e => (e.status, e.endDt, e.rc, e.executionHash, e.dockerImageHash, e.resultsClonedFrom)

    val maybeDate = if (callStatus.isStarting || callStatus.isTerminal) Option(new Date().toTimestamp) else None

    // If this call represents a call caching hit, find the execution ID for the call from which results were cloned and
    // wrap that in an `Option`.
    // If this wasn't a call caching hit just return `DBIO.successful(None)`, `None` lifted into `DBIO`.
    val findResultsClonedFromId = callStatus.resultsClonedFrom map { backendCall =>
      for {
        workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(backendCall.workflowDescriptor.id.toString).result.head
        execution <- dataAccess.executionsByWorkflowExecutionIdAndCallFqnAndIndex(
          workflowExecutionResult.workflowExecutionId.get, backendCall.key.scope.fullyQualifiedName, backendCall.key.index.fromIndex).result.head
      } yield Option(execution.executionId.get)
    } getOrElse DBIO.successful(None)

    for {
      workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.head
      executions = dataAccess.executionsByWorkflowExecutionIdAndScopeKeys(workflowExecutionResult.workflowExecutionId.get, scopeKeys)
      clonedFromId <- findResultsClonedFromId
      overallHash = callStatus.hash map { _.overallHash }
      dockerHash = callStatus.hash flatMap { _.dockerHash }
      count <- executions.map(projectionFn).update((callStatus.executionStatus.toString, maybeDate, callStatus.returnCode, overallHash, dockerHash, clonedFromId))
      scopeSize = scopeKeys.size
      _ = require(count == scopeSize, s"Execution update count $count did not match scopes size $scopeSize")
    } yield ()
  }

  override def setStatus(workflowId: WorkflowId, scopeKeys: Traversable[ExecutionDatabaseKey],
                         callStatus: CallStatus): Future[Unit] = {
    if (scopeKeys.isEmpty) Future.successful(()) else runTransaction(setStatusAction(workflowId, scopeKeys, callStatus))
  }

  override def getExecutions(id: WorkflowId): Future[Traversable[Execution]] = {
    val action = dataAccess.executionsByWorkflowExecutionUuid(id.toString).result

    runTransaction(action)
  }

  override def getExecutionsForRestart(id: WorkflowId): Future[Traversable[Execution]] = {
    val action = dataAccess.executionsForRestartByWorkflowExecutionUuid(id.toString).result

    runTransaction(action)
  }

  override def getExecutionsWithResuableResultsByHash(hash: String): Future[Traversable[Execution]] = {
    val action = dataAccess.executionsWithReusableResultsByExecutionHash(hash).result

    runTransaction(action)
  }

  override def getWorkflowExecution(workflowId: WorkflowId): Future[WorkflowExecution] = {
    val action = dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.headOption

    runTransaction(action) map { _.getOrElse(throw new NoSuchElementException(s"Workflow $workflowId not found.")) }
  }

  override def getWorkflowExecutionAux(id: WorkflowId): Future[WorkflowExecutionAux] = {
    val action = dataAccess.workflowExecutionAuxesByWorkflowExecutionUuid(id.toString).result.headOption

    runTransaction(action) map { _.getOrElse(throw new NoSuchElementException(s"No workflow execution aux found for ID '$id'.")) }
  }

  override def getAllInputs(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]] = {
    val action = dataAccess.symbolsByWorkflowExecutionUuidAndIo(workflowId.toString, IoInput).result

    runTransaction(action) map toSymbolStoreEntries
  }

  override def getAllOutputs(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]] = {
    val action = dataAccess.symbolsByWorkflowExecutionUuidAndIo(workflowId.toString, IoOutput).result

    runTransaction(action) map toSymbolStoreEntries
  }

  override def jesJobInfo(id: WorkflowId): Future[Map[ExecutionDatabaseKey, JesJob]] = {
    val action = for {
      executionAndJob <- dataAccess.jesJobsWithExecutionsByWorkflowExecutionUuid(id.toString).result
    } yield executionAndJob

    runTransaction(action) map { results =>
      results map { case (execution, job) =>
          ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex) -> job
      } toMap
    }
  }

  override def localJobInfo(id: WorkflowId): Future[Map[ExecutionDatabaseKey, LocalJob]] = {
    val action = for {
      executionAndJob <- dataAccess.localJobsWithExecutionsByWorkflowExecutionUuid(id.toString).result
    } yield executionAndJob

    runTransaction(action) map { results =>
      results map { case (execution, job) =>
        ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex) -> job
      } toMap
    }
  }

  override def sgeJobInfo(id: WorkflowId): Future[Map[ExecutionDatabaseKey, SgeJob]] = {
    val action = for {
      executionAndJob <- dataAccess.sgeJobsWithExecutionsByWorkflowExecutionUuid(id.toString).result
    } yield executionAndJob

    runTransaction(action) map { results =>
      results map { case (execution, job) =>
        ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex) -> job
      } toMap
    }
  }

  override def updateWorkflowOptions(workflowId: WorkflowId, workflowOptionsJson: String): Future[Unit] = {
    val action = for {
      workflowExecution <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.id.toString).result.head
      count <- dataAccess.workflowOptionsFromWorkflowId(workflowExecution.workflowExecutionId.get).update(workflowOptionsJson.toClob)
      _ = require(count == 1, s"Unexpected workflow aux update count $count")
    } yield ()

    runTransaction(action)
  }

  override def resetNonResumableJesExecutions(workflowId: WorkflowId): Future[Unit] = {
    // These executions have no corresponding recorded operation ID and are therefore not resumable.
    def collectNonResumableDatabaseKeys(executionsAndJobs: Seq[(Execution, JesJob)]): Seq[ExecutionDatabaseKey] = {
      executionsAndJobs collect {
        case (execution, job) if execution.status.toExecutionStatus == ExecutionStatus.Running && job.jesRunId.isEmpty =>
          ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex)
      }
    }

    val action = for {
      executionsAndJobs <- dataAccess.jesJobsWithExecutionsByWorkflowExecutionUuid(workflowId.toString).result
      nonResumableDatabaseKeys = collectNonResumableDatabaseKeys(executionsAndJobs)
      _ <- setStatusAction(workflowId, nonResumableDatabaseKeys, CallStatus(ExecutionStatus.NotStarted, None, None, None))
    } yield ()

    runTransaction(action)
  }

  override def findResumableJesExecutions(workflowId: WorkflowId): Future[Map[ExecutionDatabaseKey, JesJobKey]] = {
    // These executions have a corresponding recorded operation ID and should therefore be resumable.
    def collectResumableKeyPairs(executionsAndJobs: Traversable[(Execution, JesJob)]): Traversable[(ExecutionDatabaseKey, JesJobKey)] = {
      executionsAndJobs collect {
        case (execution, job) if execution.status.toExecutionStatus == ExecutionStatus.Running && job.jesRunId.nonEmpty =>
          (ExecutionDatabaseKey(execution.callFqn, execution.index.toIndex), JesJobKey(job.jesRunId.get))
      }
    }

    val action = for {
      executionsAndJobs <- dataAccess.jesJobsWithExecutionsByWorkflowExecutionUuid(workflowId.toString).result
      resumableKeyPairs = collectResumableKeyPairs(executionsAndJobs)
    } yield resumableKeyPairs

    runTransaction(action) map { _.toMap }
  }

  override def queryWorkflows(queryParameters: WorkflowQueryParameters): Future[WorkflowQueryResponse] = {
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

  override def updateCallCaching(parameters: CallCachingParameters): Future[Int] = {
    // Figure out which of the three possible queries to use based on whether a call has been specified and
    // if so whether an index has been specified.
    val executionQuery: (Int) => Query[dataAccess.Executions, Execution, Seq] = {
      (parameters.callKey, parameters.callKey flatMap { _.index }) match {
        case (Some(key), Some(idx)) => dataAccess.executionsByWorkflowExecutionIdAndCallFqnAndIndex(_: Int, key.fqn, idx).extract
        case (Some(key), None) => dataAccess.executionsByWorkflowExecutionIdAndCallFqn(_: Int, key.fqn).extract
        case _ => dataAccess.executionsByWorkflowExecutionId(_: Int).extract
      }
    }

    val action = for {
      workflowExecution <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(parameters.workflowId.id.toString).result.head
      count <- executionQuery(workflowExecution.workflowExecutionId.get).map(_.allowsResultReuse).update(parameters.allow)
    } yield count

    runTransaction(action)
  }
}
