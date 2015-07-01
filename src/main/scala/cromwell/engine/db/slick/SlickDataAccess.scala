package cromwell.engine.db.slick

import java.sql.Timestamp
import java.util.{Date, UUID}

import _root_.slick.util.ConfigExtensionMethods._
import com.typesafe.config.Config
import cromwell.binding._
import cromwell.binding.types.WdlType
import cromwell.binding.values.WdlValue
import cromwell.engine._
import cromwell.engine.backend.Backend
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.db._

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
}

// TODO: Still needs more refactoring

// Example: Instead of nested .filter(), breakup and compose actions, then use for-comprehension "if" style
// - Composing: https://github.com/slick/slick/blob/master/slick-testkit/src/main/scala/com/typesafe/slick/testkit/tests/TransactionTest.scala
// - .filter(): http://slick.typesafe.com/doc/3.0.0/queries.html#sorting-and-filtering
// - for based if: http://slick.typesafe.com/doc/3.0.0/queries.html#monadic-joins
// While refactoring, ran into an error: "polymorphic expression cannot be instantiated to expected type"
// Based on google, may be a problem with implicits, requiring explicit specification of some types.

// Still needs compiled slick-queries too http://slick.typesafe.com/doc/3.0.0/queries.html#compiled-queries

class SlickDataAccess(databaseConfig: Config, val dataAccess: DataAccessComponent) extends DataAccess {

  def this(databaseConfig: Config) = this(
    databaseConfig,
    new DataAccessComponent(databaseConfig.getString("slick.driver")))

  def this() = this(
    DatabaseConfig.databaseConfig)

  // NOTE: Used for slick flatMap. May switch to custom ExecutionContext the future
  private implicit val executionContext = ExecutionContext.global

  import SlickDataAccess._

  // Allows creation of a Database, plus implicits for running transactions
  import dataAccess.driver.api._

  // NOTE: if you want to refactor database is inner-class type: this.dataAccess.driver.backend.DatabaseFactory
  val database = Database.forConfig("", databaseConfig)

  // Possibly create the database
  {
    if (databaseConfig.getBooleanOr("slick.createSchema")) {
      Await.result(database.run(dataAccess.schema.create), Duration.Inf)
    }
  }

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
        workflowInfo.wdlSource,
        workflowInfo.wdlJson)

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
        symbol.wdlValue.map(_.toRawString.get))
    }
  }

  // Converts the Traversable[Call] to Seq[DBIOAction[]] that insert the correct rows
  private def toCallActions(workflowExecution: WorkflowExecution, backend: Backend,
                            calls: Traversable[Call]): Seq[DBIO[Unit]] = {
    def toWorfklowExecutionCallAction(call: Call) = toCallAction(workflowExecution, backend, call)
    calls.toSeq map toWorfklowExecutionCallAction
  }

  // Converts a single Call to a composite DBIOAction[] that inserts the correct rows
  private def toCallAction(workflowExecution: WorkflowExecution, backend: Backend,
                           call: Call): DBIO[Unit] = {
    for {
    // Insert an execution row
      executionInsert <- dataAccess.executionsAutoInc +=
        new Execution(
          workflowExecution.workflowExecutionId.get,
          call.taskFqn,
          ExecutionStatus.NotStarted.toString)

      // Depending on the backend, insert a job specific row
      _ <- backend match {
        case _: LocalBackend =>
          dataAccess.localJobsAutoInc +=
            new LocalJob(
              executionInsert.executionId.get,
              None,
              None)
        case null =>
          throw new NotImplementedError("Backend is null")
        case unknown =>
          throw new NotImplementedError("Unknown backend: " + backend.getClass)
      }
    } yield ()
  }

  override def getWorkflowState(workflowId: WorkflowId): Future[Option[WorkflowState]] = {

    val action = for {

      workflowExecutionResultOption <- dataAccess.workflowExecutions.filter(
        _.workflowExecutionUuid === workflowId.toString).result.headOption

      workflowState = workflowExecutionResultOption map { workflowExecutionResult =>
        WorkflowState.fromString(workflowExecutionResult.status)
      }

    } yield workflowState

    runTransaction(action)
  }

  override def getExecutionStatuses(workflowId: WorkflowId): Future[Map[FullyQualifiedName, CallStatus]] = {

    val action = for {

      workflowExecutionResult <- dataAccess.workflowExecutions.filter(
        _.workflowExecutionUuid === workflowId.toString).result.head

      executionResults <- dataAccess.executions.filter(symbol =>
        symbol.workflowExecutionId === workflowExecutionResult.workflowExecutionId).result

      executionStatuses = executionResults.map(executionResult =>
        executionResult.callFqn -> ExecutionStatus.withName(executionResult.status)
      ).toMap

    } yield executionStatuses

    runTransaction(action)
  }

  override def getWorkflowsByState(states: Traversable[WorkflowState]): Future[Traversable[WorkflowInfo]] = {

    val action = for {

      workflowExecutionResults <- dataAccess.workflowExecutions.filter(
        _.status inSet states.map(_.toString)).result

      workflowInfos <- DBIO.sequence(
        workflowExecutionResults map { workflowExecutionResult =>
          val workflowExecutionAuxResult = dataAccess.workflowExecutionAuxes.filter(
            _.workflowExecutionId === workflowExecutionResult.workflowExecutionId).result.head

          workflowExecutionAuxResult map { workflowExecutionAux =>
            new WorkflowInfo(
              UUID.fromString(workflowExecutionResult.workflowExecutionUuid),
              workflowExecutionAux.wdlSource,
              workflowExecutionAux.jsonInputs)
          }
        }
      )

    } yield workflowInfos

    runTransaction(action)
  }

  override def updateWorkflowState(workflowId: WorkflowId, workflowState: WorkflowState): Future[Unit] = {
    val query = for {
      workflowExecution <- dataAccess.workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowId.toString
    } yield workflowExecution

    val action = for {
      count <- query.map(_.status).update(workflowState.toString)
      _ = require(count == 1, s"Unexpected workflow execution update count $count")
    } yield ()

    runTransaction(action)
  }

  override def getExecutionBackendInfo(workflowId: WorkflowId, call: Call): Future[CallBackendInfo] = {

    val action = for {

      workflowExecutionResult <- dataAccess.workflowExecutions.filter(
        _.workflowExecutionUuid === workflowId.toString).result.head

      executionResult <- dataAccess.executions.filter(execution =>
        execution.workflowExecutionId === workflowExecutionResult.workflowExecutionId &&
          execution.callFqn === call.fullyQualifiedName).result.head

      localJobResultOption <- dataAccess.localJobs.filter(localJob =>
        localJob.executionId === executionResult.executionId).result.headOption

      jesJobResultOption <- dataAccess.jesJobs.filter(jesJob =>
        jesJob.executionId === executionResult.executionId).result.headOption

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
          throw new NotImplementedError(
            s"Unknown backend from db for (uuid, fqn): " +
              s"($workflowId, ${call.fullyQualifiedName})")
      }

    } yield backendInfo

    runTransaction(action)
  }

  override def updateExecutionBackendInfo(workflowId: WorkflowId,
                                          call: Call,
                                          backendInfo: CallBackendInfo): Future[Unit] = {
    val action = for {

      workflowExecutionResult <- dataAccess.workflowExecutions.filter(
        _.workflowExecutionUuid === workflowId.toString).result.head

      executionQuery = dataAccess.executions.filter(execution =>
        execution.workflowExecutionId === workflowExecutionResult.workflowExecutionId &&
          execution.callFqn === call.fullyQualifiedName)

      executionUpdate <- executionQuery.map(_.status).update(backendInfo.status.toString)

      _ = require(executionUpdate == 1, s"Unexpected execution update count $executionUpdate")

      executionResult <- executionQuery.result.head

      backendUpdate <- backendInfo match {
        case localBackendInfo: LocalCallBackendInfo =>
          dataAccess.localJobs.filter(
            _.executionId === executionResult.executionId
          ).map(cols => (
            cols.pid,
            cols.rc
            )).update(
              localBackendInfo.processId,
              localBackendInfo.resultCode
            )

        case jesBackendInfo: JesCallBackendInfo =>
          dataAccess.jesJobs.filter(
            _.executionId === executionResult.executionId
          ).map(cols => (
            cols.jesId,
            cols.jesStatus
            )).update(
              jesBackendInfo.jesId,
              jesBackendInfo.jesStatus
            )

        case null =>
          throw new NotImplementedError("Backend info is null")
      }

      _ = require(backendUpdate == 1, s"Unexpected backend update count $backendUpdate")

    } yield ()

    runTransaction(action)
  }

  /** Get all inputs for the scope of this call. */
  override def getInputs(workflowId: WorkflowId, call: Call): Future[Traversable[SymbolStoreEntry]] = {
    require(call != null, "call cannot be null")
    getSymbols(workflowId, IoInput, Option(call))
  }

  /** Get all outputs for the scope of this call. */
  override def getOutputs(workflowId: WorkflowId, call: Call): Future[Traversable[SymbolStoreEntry]] = {
    require(call != null, "call cannot be null")
    getSymbols(workflowId, IoOutput, Option(call))
  }

  /** Returns all outputs for this workflowId */
  override def getOutputs(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]] = {
    getSymbols(workflowId, IoOutput, None)
  }

  private def getSymbols(workflowId: WorkflowId, ioValue: IoValue,
                         callOption: Option[Call] = None): Future[Traversable[SymbolStoreEntry]] = {

    val action = for {

      workflowExecutionResult <- dataAccess.workflowExecutions.filter(
        _.workflowExecutionUuid === workflowId.toString).result.head

      symbolResults <- dataAccess.symbols.filter(symbol =>
        symbol.workflowExecutionId === workflowExecutionResult.workflowExecutionId &&
          symbol.io === ioValue && {
          callOption match {
            case Some(call) => symbol.scope === call.fullyQualifiedName
            case None => true: Rep[Boolean]
          }
        }).result

      symbolStoreEntries = symbolResults map { symbolResult =>
        val wdlType = WdlType.fromWdlString(symbolResult.wdlType)
        val symbolStoreKey = new SymbolStoreKey(
          symbolResult.scope,
          symbolResult.name,
          Option(symbolResult.iteration).filterNot(_ == IterationNone),
          input = symbolResult.io == IoInput) // input = true, if db contains "INPUT"
        new SymbolStoreEntry(
          symbolStoreKey,
          wdlType,
          symbolResult.wdlValue.map(wdlType.coerceRawValue(_).get))

      }

    } yield symbolStoreEntries

    runTransaction(action)
  }

  /** Should fail if a value is already set.  The keys in the Map are locally qualified names. */
  override def setOutputs(workflowId: WorkflowId, call: Call, callOutputs: Map[String, WdlValue]): Future[Unit] = {

    val query = for {
      workflowExecution <- dataAccess.workflowExecutions
      if workflowExecution.workflowExecutionUuid === workflowId.toString
    } yield workflowExecution

    val action = for {
      workflowExecution <- query.result.head
      _ <- dataAccess.symbolsAutoInc ++= callOutputs map {
        case (symbolLocallyQualifiedName, wdlValue) =>
          new Symbol(
            workflowExecution.workflowExecutionId.get,
            call.fullyQualifiedName,
            call.fullyQualifiedName + "." + symbolLocallyQualifiedName,
            IterationNone,
            IoOutput,
            wdlValue.wdlType.toWdlString,
            Option(wdlValue.toRawString.get))
      }
    } yield ()

    runTransaction(action)
  }

  override def setStatus(workflowId: WorkflowId, calls: Traversable[Call], callStatus: CallStatus): Future[Unit] = {

    val action = for {
      workflowExecutionResult <- dataAccess.workflowExecutions.filter(
        _.workflowExecutionUuid === workflowId.toString).result.head

      count <- dataAccess.executions.filter(execution =>
        (execution.workflowExecutionId === workflowExecutionResult.workflowExecutionId.get) &&
          (execution.callFqn inSet calls.map(_.fullyQualifiedName))
      ).map(_.status).update(callStatus.toString)

      callSize = calls.size
      _ = require(count == callSize, s"Execution update count $count did not match calls size $callSize")

    } yield ()

    runTransaction(action)
  }
}
