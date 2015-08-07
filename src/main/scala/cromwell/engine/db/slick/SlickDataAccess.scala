package cromwell.engine.db.slick

import java.sql.{Clob, Timestamp}
import java.util.{Date, UUID}
import javax.sql.rowset.serial.SerialClob

import _root_.slick.util.ConfigExtensionMethods._
import com.typesafe.config.{Config, ConfigFactory, ConfigValueFactory}
import cromwell.binding._
import cromwell.binding.types.WdlType
import cromwell.engine._
import cromwell.engine.backend.Backend
import cromwell.engine.backend.jes.JesBackend
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.db._
import org.slf4j.LoggerFactory

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, ExecutionContext, Future}
import scala.language.implicitConversions

object SlickDataAccess {
  type IoValue = String
  val IoInput = "INPUT"
  val IoOutput = "OUTPUT"
  val IterationNone = -1 // "It's a feature" https://bugs.mysql.com/bug.php?id=8173

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
    def urlKey = if (config.hasPath("url")) "url" else "properties.url"

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

class SlickDataAccess(databaseConfig: Config, val dataAccess: DataAccessComponent) extends DataAccess {

  def this(databaseConfig: Config) = this(
    databaseConfig,
    new DataAccessComponent(databaseConfig.getString("slick.driver")))

  def this() = this(DatabaseConfig.databaseConfig)

  // NOTE: Used for slick flatMap. May switch to custom ExecutionContext the future
  private implicit val executionContext = ExecutionContext.global

  import SlickDataAccess._

  // Allows creation of a Database, plus implicits for running transactions
  import dataAccess.driver.api._

  // NOTE: if you want to refactor database is inner-class type: this.dataAccess.driver.backend.DatabaseFactory
  private val configWithUniqueSchema = databaseConfig.withUniqueSchema
  val database = Database.forConfig("", configWithUniqueSchema)

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
    if (databaseConfig.getBooleanOr("slick.createSchema")) {
      Await.result(database.run(dataAccess.schema.create), Duration.Inf)
    }
  }

  override def shutdown() = database.shutdown

  // Run action with an outer transaction
  private def runTransaction[R](action: DBIOAction[R, _ <: NoStream, _ <: Effect]): Future[R] = {
    database.run(action.transactionally)
  }

  /**
   * Creates a row in each of the backend-info specific tables for each call in `calls` corresponding to the backend
   * `backend`.  Or perhaps defer this?
   */
  override def createWorkflow(workflowInfo: WorkflowInfo,
                              workflowInputs: Traversable[SymbolStoreEntry],
                              calls: Traversable[Call],
                              backend: Backend): Future[Unit] = {

    val action = for {

      workflowExecutionInsert <- dataAccess.workflowExecutionsAutoInc +=
        new WorkflowExecution(
          workflowInfo.workflowId.toString,
          WorkflowSubmitted.toString,
          new Date().toTimestamp)

      _ <- dataAccess.workflowExecutionAuxesAutoInc += new WorkflowExecutionAux(
        workflowExecutionInsert.workflowExecutionId.get,
        workflowInfo.wdlSource.toClob,
        workflowInfo.wdlJson.toClob)

      symbolInsert <- dataAccess.symbolsAutoInc ++= toSymbols(workflowExecutionInsert, workflowInputs)

      // NOTE: Don't use DBIO.seq for **transforming** sequences
      // - DBIO.seq(mySeq: _*) runs *any* items in sequence, but converts Seq[ DBIOAction[_] ] to DBIOAction[ Unit ]
      // - DBIO.sequence(mySeq) converts Seq[ DBIOAction[R] ] to DBIOAction[ Seq[R] ]
      // - DBIO.fold(mySeq, init) converts Seq[ DBIOAction[R] ] to DBIOAction[R]

      _ <- DBIO.sequence(toCallActions(workflowExecutionInsert, backend, calls))

    } yield ()

    runTransaction(action)
  }

  // Converts the SymbolStoreEntry to Symbols. Does not create the action to do the insert.
  private def toSymbols(workflowExecution: WorkflowExecution,
                        symbolStoreEntries: Traversable[SymbolStoreEntry]): Seq[Symbol] = {
    symbolStoreEntries.toSeq map { symbol =>
      new Symbol(
        workflowExecution.workflowExecutionId.get,
        symbol.key.scope,
        symbol.key.name,
        symbol.key.iteration.getOrElse(IterationNone),
        if (symbol.isInput) IoInput else IoOutput,
        symbol.wdlType.toWdlString,
        symbol.wdlValue.map(_.toWdlString.toClob))
    }
  }

  // Converts the Traversable[Call] to Seq[DBIOAction[]] that insert the correct rows
  private def toCallActions(workflowExecution: WorkflowExecution, backend: Backend,
                            calls: Traversable[Call]): Seq[DBIO[Unit]] = {
    def toWorkflowExecutionCallAction(call: Call) = toCallAction(workflowExecution, backend, call)
    calls.toSeq map toWorkflowExecutionCallAction
  }

  // Converts a single Call to a composite DBIOAction[] that inserts the correct rows
  private def toCallAction(workflowExecution: WorkflowExecution, backend: Backend,
                           call: Call): DBIO[Unit] = {
    for {
    // Insert an execution row
      executionInsert <- dataAccess.executionsAutoInc +=
        new Execution(
          workflowExecution.workflowExecutionId.get,
          call.fullyQualifiedName,
          ExecutionStatus.NotStarted.toString)

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
          dataAccess.jesJobsAutoInc += new JesJob(executionInsert.executionId.get, 0, "", None)
        case null =>
          throw new IllegalArgumentException("Backend is null")
        case unknown =>
          throw new IllegalArgumentException("Unknown backend: " + backend.getClass)
      }

    } yield ()
  }

  override def getWorkflowState(workflowId: WorkflowId): Future[Option[WorkflowState]] = {

    val action = for {
      workflowExecutionStatusOption <- dataAccess.workflowExecutionStatusesByWorkflowExecutionUuid(
        workflowId.toString).result.headOption
      workflowState = workflowExecutionStatusOption map WorkflowState.fromString
    } yield workflowState

    runTransaction(action)
  }

  override def getExecutionStatuses(workflowId: WorkflowId): Future[Map[FullyQualifiedName, CallStatus]] = {

    val action = for {

    // NOTE: For now, intentionally causes query to error out instead of returning an Map.empty
      workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(
        workflowId.toString).result.head

      // Alternatively, could use a dataAccess.executionCallFqnsAndStatusesByWorkflowExecutionUuid
      executionCallFqnAndStatusResults <- dataAccess.executionCallFqnsAndStatusesByWorkflowExecutionId(
        workflowExecutionResult.workflowExecutionId.get).result

      executionStatuses = executionCallFqnAndStatusResults.toMap mapValues ExecutionStatus.withName

    } yield executionStatuses

    runTransaction(action)
  }

  override def getExecutionStatus(workflowId: WorkflowId, callFqn: FullyQualifiedName): Future[Option[CallStatus]] = {
    val action = for {
      workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(
        workflowId.toString).result.head

      executionStatuses <- dataAccess.executionStatusByWorkflowExecutionIdAndCallFqn(
        (workflowExecutionResult.workflowExecutionId.get, callFqn)).result

      maybeStatus = executionStatuses.headOption map ExecutionStatus.withName
    } yield maybeStatus
    runTransaction(action)
  }

  override def getWorkflowsByState(states: Traversable[WorkflowState]): Future[Traversable[WorkflowInfo]] = {

    val action = for {

      workflowExecutionResults <- dataAccess.workflowExecutionsByStatuses(states.map(_.toString)).result

      workflowInfos <- DBIO.sequence(
        workflowExecutionResults map { workflowExecutionResult =>

          val workflowExecutionAuxResult = dataAccess.workflowExecutionAuxesByWorkflowExecutionId(
            workflowExecutionResult.workflowExecutionId.get).result.head

          workflowExecutionAuxResult map { workflowExecutionAux =>
            new WorkflowInfo(
              UUID.fromString(workflowExecutionResult.workflowExecutionUuid),
              workflowExecutionAux.wdlSource.toRawString,
              workflowExecutionAux.jsonInputs.toRawString)
          }
        }
      )

    } yield workflowInfos

    runTransaction(action)
  }

  override def updateWorkflowState(workflowId: WorkflowId, workflowState: WorkflowState): Future[Unit] = {
    val action = for {
      count <- dataAccess.workflowExecutionStatusesByWorkflowExecutionUuid(
        workflowId.toString).update(workflowState.toString)
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

      jobResultOption = localJobResultOption orElse jesJobResultOption
      backendInfo = jobResultOption match {
        case Some(localJobResult: LocalJob) =>
          new LocalCallBackendInfo(
            ExecutionStatus.withName(executionResult.status),
            localJobResult.pid,
            localJobResult.rc)
        case Some(jesJobResult: JesJob) =>
          new JesCallBackendInfo(
            ExecutionStatus.withName(executionResult.status),
            jesJobResult.jesId,
            jesJobResult.jesStatus)
        case _ =>
          throw new IllegalArgumentException(
            s"Unknown backend from db for (uuid, fqn): " +
              s"($workflowId, ${call.fullyQualifiedName})")
      }

    } yield backendInfo

    runTransaction(action)
  }

  override def updateExecutionBackendInfo(workflowId: WorkflowId,
                                          call: Call,
                                          backendInfo: CallBackendInfo): Future[Unit] = {

    require(backendInfo != null, "backend info is null")

    val action = for {

      executionResult <- dataAccess.executionsByWorkflowExecutionUuidAndCallFqn(
        workflowId.toString, call.fullyQualifiedName).result.head

      executionStatusQuery = dataAccess.executionStatusesByExecutionId(
        executionResult.executionId.get)

      executionUpdate <- executionStatusQuery.update(
        backendInfo.status.toString)

      _ = require(executionUpdate == 1, s"Unexpected execution update count $executionUpdate")

      backendUpdate <- backendInfo match {
        case localBackendInfo: LocalCallBackendInfo =>
          dataAccess.localJobPidsAndRcsByExecutionId(
            executionResult.executionId.get).update(
              localBackendInfo.processId,
              localBackendInfo.resultCode)

        case jesBackendInfo: JesCallBackendInfo =>
          dataAccess.jesJobIdsAndJesStatusesByExecutionId(
            executionResult.executionId.get).update(
              jesBackendInfo.jesId,
              jesBackendInfo.jesStatus
            )
      }

      _ = require(backendUpdate == 1, s"Unexpected backend update count $backendUpdate")

    } yield ()

    runTransaction(action)
  }

  private def toSymbolStoreEntry(symbolResult: Symbol) = {
    val wdlType = WdlType.fromWdlString(symbolResult.wdlType)
    new SymbolStoreEntry(
      new SymbolStoreKey(
        symbolResult.scope,
        symbolResult.name,
        Option(symbolResult.iteration).filterNot(_ == IterationNone),
        input = symbolResult.io == IoInput // input = true, if db contains "INPUT"
      ),
      wdlType,
      symbolResult.wdlValue map {value => wdlType.fromWdlString(value.toRawString)}
    )
  }

  override def getFullyQualifiedName(workflowId: WorkflowId, fqn: FullyQualifiedName): Future[Traversable[SymbolStoreEntry]] = {
    val Array(scope, varName) = fqn.split("\\.(?=[^\\.]+$)") // e.g. "a.b.c.d" => Seq("a.b.c", "d")
    val action = for {
      symbolResults <- dataAccess.symbolsByScopeAndName(workflowId.toString, scope, varName).result
      symbolStoreEntries = symbolResults map toSymbolStoreEntry
    } yield symbolStoreEntries

    runTransaction(action)
  }

  override def getAll(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]] = {
    val x = dataAccess.allSymbols(workflowId.toString).result
    val action = for {
      symbolResults <- x
      symbolStoreEntries = symbolResults map toSymbolStoreEntry
    } yield symbolStoreEntries

    runTransaction(action)
  }

  /** Get all inputs for the scope of this call. */
  override def getInputs(workflowId: WorkflowId, call: Call): Future[Traversable[SymbolStoreEntry]] = {
    require(call != null, "call cannot be null")
    getSymbols(workflowId, IoInput, Option(call.fullyQualifiedName))
  }

  /** Get all outputs for the scope of this call. */
  override def getOutputs(workflowId: WorkflowId, callFqn: FullyQualifiedName): Future[Traversable[SymbolStoreEntry]] = {
    require(callFqn != null, "callFqn cannot be null")
    getSymbols(workflowId, IoOutput, Option(callFqn))
  }

  /** Returns all outputs for this workflowId */
  override def getOutputs(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]] = {
    getSymbols(workflowId, IoOutput, None)
  }

  private def getSymbols(workflowId: WorkflowId, ioValue: IoValue,
                         callFqnOption: Option[FullyQualifiedName] = None): Future[Traversable[SymbolStoreEntry]] = {

    val action = for {
      symbolResults <- dataAccess.symbolsByWorkflowExecutionUuidAndIoAndMaybeScope(
        workflowId.toString, ioValue, callFqnOption
      ).result
      symbolStoreEntries = symbolResults map toSymbolStoreEntry
    } yield symbolStoreEntries

    runTransaction(action)
  }

  /** Should fail if a value is already set.  The keys in the Map are locally qualified names. */
  override def setOutputs(workflowId: WorkflowId, call: Call, callOutputs: WorkflowOutputs): Future[Unit] = {

    val action = for {
      workflowExecution <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.head
      _ <- dataAccess.symbolsAutoInc ++= callOutputs map {
        case (symbolLocallyQualifiedName, wdlValue) =>
          new Symbol(
            workflowExecution.workflowExecutionId.get,
            call.fullyQualifiedName,
            symbolLocallyQualifiedName,
            IterationNone,
            IoOutput,
            wdlValue.wdlType.toWdlString,
            Option(wdlValue.toWdlString.toClob))
      }
    } yield ()

    runTransaction(action)
  }

  override def setStatus(workflowId: WorkflowId, callFqns: Traversable[FullyQualifiedName],
                         callStatus: CallStatus): Future[Unit] = {

    val action = for {
      workflowExecutionResult <- dataAccess.workflowExecutionsByWorkflowExecutionUuid(workflowId.toString).result.head

      count <- dataAccess.executionStatusesByWorkflowExecutionIdAndCallFqns(
        workflowExecutionResult.workflowExecutionId.get, callFqns).update(callStatus.toString)

      callSize = callFqns.size
      _ = require(count == callSize, s"Execution update count $count did not match calls size $callSize")

    } yield ()

    runTransaction(action)
  }
}
