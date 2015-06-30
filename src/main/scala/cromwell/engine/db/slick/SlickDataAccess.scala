package cromwell.engine.db.slick

import java.util.{Date, UUID}

import com.typesafe.config.Config
import cromwell.binding._
import cromwell.binding.types.WdlType
import cromwell.binding.values.{WdlPrimitive, WdlValue}
import cromwell.engine._
import cromwell.engine.backend.Backend
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.db._

import scala.concurrent.{ExecutionContext, Future}
import scala.language.implicitConversions

object SlickDataAccess {
  implicit class DateToTimestamp(val date: java.util.Date) extends AnyVal {
    def toTimestamp = new java.sql.Timestamp(date.getTime)
  }
}

// TODO: Still needs more refactoring

// Example: Instead of nested .filter(), breakup and compose actions, then use for-comprehension "if" style
// - Composing: https://github.com/slick/slick/blob/master/slick-testkit/src/main/scala/com/typesafe/slick/testkit/tests/TransactionTest.scala
// - .filter(): http://slick.typesafe.com/doc/3.0.0/queries.html#sorting-and-filtering
// - for based if: http://slick.typesafe.com/doc/3.0.0/queries.html#monadic-joins
// While refactoring, ran into an error: "polymorphic expression cannot be instantiated to expected type"
// Based on google, may be a problem with implicits, requiring explicit specification of some types.

class SlickDataAccess(databaseConfig: Config) extends DataAccess {

  def this() = this(DatabaseConfig.databaseConfig)

  private val IoInput = "INPUT"
  private val IoOutput = "OUTPUT"

  import SlickDataAccess._

  val dataAccess = new DataAccessComponent(databaseConfig.getString("slick.driver"))

  import dataAccess.driver.api._

  // TODO: Used for flatMap. Should we build another ec instead of using global?
  private implicit val executionContext = ExecutionContext.global

  // WdlValues have multiple ways they coerce back. This currently decoerces for a proper round trip.
  // TODO: Refactor into the type?
  private def wdlValueToString(wdlValue: WdlValue) = {
    wdlValue match {
      case wdlPrimitive: WdlPrimitive => wdlPrimitive.asString
      case other => other.toString
    }
  }

  val database = Database.forConfig("", databaseConfig)

  // Check the database connection. Can be run before operations that actually use the database.
  def isValidConnection(timeoutSecs: Int): Future[Boolean] = {
    database.run(SimpleDBIO(_.connection.isValid(timeoutSecs))) recover { case _ => false }
  }

  // Lazily, possibly create the database, returning the result in the future
  private lazy val createDatabase: Future[Unit] = {
    import _root_.slick.util.ConfigExtensionMethods._
    if (databaseConfig.getBooleanOr("slick.createSchema")) {
      database.run(dataAccess.schema.create)
    } else Future.successful(())
  }

  // Run action with an outer transaction. Also if we need to, create the in memory database!
  private def runTransaction[R](action: DBIOAction[R, _ <: NoStream, _ <: Effect]): Future[R] = {
    for {
      _ <- createDatabase // If not created
      result <- database.run(action.transactionally)
    } yield result
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

      symbolInsert <- dataAccess.symbolsAutoInc ++=
        workflowInputs.toSeq map { symbol =>
          new Symbol(
            workflowExecutionInsert.workflowExecutionId.get,
            symbol.key.scope,
            symbol.key.name,
            symbol.key.iteration,
            if (symbol.isInput) IoInput else IoOutput,
            symbol.wdlType.toWdlString,
            symbol.wdlValue.map(wdlValueToString).orNull)
        }

      // NOTE: Don't use DBIO.seq for **transforming** sequences
      // - DBIO.seq(mySeq: _*) runs *any* items in sequence, but converts Seq[ DBIOAction[_] ] to DBIOAction[ Unit ]
      // - DBIO.sequence(mySeq) converts Seq[ DBIOAction[R] ] to DBIOAction[ Seq[R] ]
      // - DBIO.fold(mySeq, init) converts Seq[ DBIOAction[R] ] to DBIOAction[R]

      _ <- DBIO.sequence(calls.toSeq map {
        call => for {
          executionInsert <- dataAccess.executionsAutoInc +=
            new Execution(
              workflowExecutionInsert.workflowExecutionId.get,
              call.taskFqn,
              ExecutionStatus.NotStarted.toString)

          _ <- backend match {
            case _: LocalBackend =>
              dataAccess.localJobsAutoInc +=
                new LocalJob(
                  executionInsert.executionId.get,
                  -1,
                  -1)
            case null =>
              throw new NotImplementedError("Backend is null")
            case unknown =>
              throw new NotImplementedError("Unknown backend: " + backend.getClass)
          }
        } yield ()
      })

    } yield ()

    runTransaction(action)
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

      backendInfo =
      if (localJobResultOption.isDefined) {
        new LocalCallBackendInfo(
          ExecutionStatus.withName(executionResult.status),
          localJobResultOption.get.pid,
          localJobResultOption.get.rc)
      } else if (jesJobResultOption.isDefined) {
        new JesCallBackendInfo(
          ExecutionStatus.withName(executionResult.status),
          jesJobResultOption.get.jesId,
          jesJobResultOption.get.jesStatus)
      } else {
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

    val action = for {

      workflowExecutionResult <- dataAccess.workflowExecutions.filter(
        _.workflowExecutionUuid === workflowId.toString).result.head

      symbolResults <- dataAccess.symbols.filter(symbol =>
        symbol.workflowExecutionId === workflowExecutionResult.workflowExecutionId &&
          symbol.io === IoInput &&
          symbol.scope === call.fullyQualifiedName).result

      symbolStoreEntries = symbolResults map { symbolResult =>
        val wdlType = WdlType.fromWdlString(symbolResult.wdlType)
        val symbolStoreKey = new SymbolStoreKey(
          symbolResult.scope,
          symbolResult.name,
          symbolResult.iteration,
          symbolResult.io == IoInput)
        new SymbolStoreEntry(
          symbolStoreKey,
          wdlType,
          Option(symbolResult.wdlValue).map(wdlType.coerceRawValue(_).get))

      }
    } yield symbolStoreEntries

    runTransaction(action)
  }

  /** Get all outputs for the scope of this call. */
  override def getOutputs(workflowId: WorkflowId, call: Call): Future[Traversable[SymbolStoreEntry]] = {

    val action = for {

      workflowExecutionResult <- dataAccess.workflowExecutions.filter(
        _.workflowExecutionUuid === workflowId.toString).result.head

      symbolResults <- dataAccess.symbols.filter(symbol =>
        symbol.workflowExecutionId === workflowExecutionResult.workflowExecutionId &&
          symbol.io === IoOutput &&
          symbol.scope === call.fullyQualifiedName).result

      symbolStoreEntries = symbolResults map { symbolResult =>
        val wdlType = WdlType.fromWdlString(symbolResult.wdlType)
        val symbolStoreKey = new SymbolStoreKey(
          symbolResult.scope,
          symbolResult.name,
          symbolResult.iteration,
          symbolResult.io == IoInput)
        new SymbolStoreEntry(
          symbolStoreKey,
          wdlType,
          Option(symbolResult.wdlValue).map(wdlType.coerceRawValue(_).get))

      }

    } yield symbolStoreEntries

    runTransaction(action)
  }

  /** Returns all outputs for this workflowId */
  override def getOutputs(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]] = {

    val action = for {

      workflowExecutionResult <- dataAccess.workflowExecutions.filter(
        _.workflowExecutionUuid === workflowId.toString).result.head

      symbolResults <- dataAccess.symbols.filter(symbol =>
        symbol.workflowExecutionId === workflowExecutionResult.workflowExecutionId &&
          symbol.io === IoOutput).result

      symbolStoreEntries = symbolResults map { symbolResult =>
        val wdlType = WdlType.fromWdlString(symbolResult.wdlType)
        val symbolStoreKey = new SymbolStoreKey(
          symbolResult.scope,
          symbolResult.name,
          symbolResult.iteration,
          symbolResult.io == IoInput)
        new SymbolStoreEntry(
          symbolStoreKey,
          wdlType,
          Option(symbolResult.wdlValue).map(wdlType.coerceRawValue(_).get))

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
        case (symbolFullyQualifiedName, wdlValue) =>
          new Symbol(
            workflowExecution.workflowExecutionId.get,
            call.fullyQualifiedName,
            symbolFullyQualifiedName,
            None,
            IoOutput,
            wdlValue.wdlType.toWdlString,
            wdlValueToString(wdlValue))
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
