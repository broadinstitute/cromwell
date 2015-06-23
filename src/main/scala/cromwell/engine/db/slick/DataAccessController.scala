package cromwell.engine.db.slick

import java.util.{UUID, Date}

import cromwell.binding.Call
import cromwell.binding.types.WdlType
import cromwell.binding.values.{WdlPrimitive, WdlValue}
import cromwell.engine.db._
import cromwell.engine.store.ExecutionStore.ExecutionStatus
import cromwell.engine.store.SymbolStore.SymbolStoreKey
import cromwell.engine.store.SymbolStoreEntry
import cromwell.engine.{WorkflowId, WorkflowState, WorkflowSubmitted}
import cromwell.util.DatabaseConfig

import scala.concurrent.duration.Duration
import scala.concurrent.{ExecutionContext, Await, Future}
import scala.language.implicitConversions

// TODO: Need lots of refactoring including: better futures handling, compiled queries, utility functions, etc.
object DataAccessController extends DataAccess {
  val dataAccess = new DataAccessComponent(DatabaseConfig.slickDriver)

  import dataAccess.driver.api._

  def database: Database = Database.forConfig("", DatabaseConfig.databaseConfig)

  // TODO: Pass in the execution context instead of using the global!!
  private implicit val executionContext = ExecutionContext.global

  private implicit def date2timestamp(date: java.util.Date): java.sql.Timestamp =
    new java.sql.Timestamp(date.getTime)

  // WdlValues have multiple ways they coerce back. This currently decoerces for a proper round trip.
  // TODO: Refactor into the type?
  private def wdlValueToString(wdlValue: WdlValue) = {
    wdlValue match {
      case wdlPrimitive: WdlPrimitive => wdlPrimitive.asString
      case other => other.toString
    }
  }

  // Will stamp the start_dt column with DateTime.Now.
  // Will stamp the WorkflowState as Submitted.
  override def createWorkflow(id: WorkflowId, wdlUri: String, symbols: Seq[SymbolStoreEntry]): Unit = {
    val insertFuture: Future[Unit] = database.run((for {
      workflowExecutionInsert <- dataAccess.insertWorkflowExecution(
        new WorkflowExecution(id.toString, wdlUri, WorkflowSubmitted.toString, new Date())
      )
      symbolsInserts <- dataAccess.insertSymbols(
        symbols map { symbol =>
          new Symbol(
            workflowExecutionInsert.workflowExecutionId.get,
            symbol.key.scope, symbol.key.name, symbol.key.iteration,
            if (symbol.isInput) "INPUT" else "OUTPUT",
            symbol.wdlType.toWdlString, symbol.wdlValue.map(wdlValueToString).orNull)
        }
      )
    } yield ()).transactionally)
    // TODO: Return Future[Unit] instead of Unit?
    Await.result(insertFuture, Duration.Inf)
  }

  override def updateCall(workflowId: WorkflowId, call: Call, notUsedByMethod: Option[CallStatus],
                          callInfoOption: Option[CallInfo], symbolsOption: Option[Seq[SymbolStoreEntry]]): Unit = {
    val workflowExecutionSelect = dataAccess.workflowExecutions.filter(
      _.workflowExecutionUuid === workflowId.toString).result.head

    val callInsertsOrUpdates = callInfoOption map { callInfo =>

      // Run a query to look for a call row with our FQN
      val callSelect = workflowExecutionSelect flatMap { workflowExecutionRow =>
        dataAccess.executions.filter(executionRow =>
          executionRow.workflowExecutionId === workflowExecutionRow.workflowExecutionId &&
            executionRow.callFqn === call.taskFqn
        ).result.headOption
      }

      // Check the query, and return a return inserts or updates
      callSelect flatMap {

        case Some(executionRow) =>
          // The select found an existing execution row. Do updates.
          val executionQuery = dataAccess.executions.filter(_.executionId === executionRow.executionId)
          val executionColumns = executionQuery.map(_.status)
          val executionUpdate = executionColumns.update(callInfo.status.toString)

          val callUpdate = callInfo match {
            case jesCallInfo: JesCallInfo =>
              val query = dataAccess.jesJobs.filter(_.executionId === executionRow.executionId.get)
              val columns = query.map(jobs => (jobs.jesId, jobs.jesStatus))
              columns.update(jesCallInfo.jesId, jesCallInfo.jesStatus)
            case localCallInfo: LocalCallInfo =>
              val query = dataAccess.localJobs.filter(_.executionId === executionRow.executionId.get)
              val columns = query.map(jobs => (jobs.pid, jobs.command, jobs.rc))
              columns.update(localCallInfo.processId, localCallInfo.command, localCallInfo.resultCode)
          }

          DBIO.seq(executionUpdate, callUpdate)

        case None =>
          for {
            workflowExecution <- workflowExecutionSelect
            executionInsert <- dataAccess.executionsAutoInc += new Execution(
              workflowExecution.workflowExecutionId.get, call.taskFqn, callInfo.status.toString)
            localCallInfo = callInfo.asInstanceOf[LocalCallInfo]
            callInsert <- callInfo match {
              case jesCallInfo: JesCallInfo =>
                dataAccess.jesJobsAutoInc += new JesJob(
                  executionInsert.executionId.get,
                  jesCallInfo.jesId, jesCallInfo.jesStatus)
              case localCallInfo: LocalCallInfo =>
                dataAccess.localJobsAutoInc += new LocalJob(
                  executionInsert.executionId.get,
                  localCallInfo.processId, localCallInfo.command, localCallInfo.resultCode)
            }
          } yield (executionInsert, callInsert)

      } // callSelect flatMap

    }

    val symbolInsertOrUpdates = symbolsOption map {
      _ map { symbol =>

        // Check each symbol to see if one exists with the (scope, name, iteration) combo
        val symbolSelect = workflowExecutionSelect flatMap { workflowExecutionRow =>
          dataAccess.symbols.filter(symbolRow =>
            symbolRow.workflowExecutionId === workflowExecutionRow.workflowExecutionId &&
              symbolRow.scope === symbol.key.scope &&
              symbolRow.name === symbol.key.name &&
              symbolRow.iteration === symbol.key.iteration
          ).result.headOption
        }

        symbolSelect map {
          case Some(symbolRow) =>
            // Select found a row. Update.
            val query = dataAccess.symbols.filter(_.symbolId === symbolRow.symbolId)
            val columns = query.map(_.wdlValue) // NOTE: Not allowing updating of the type, just the value
            columns.update(symbol.wdlValue.map(wdlValueToString).orNull)
          case None =>
            // Select did not find a row. Insert.
            workflowExecutionSelect flatMap { workflowExecutionRow =>
              dataAccess.symbolsAutoInc += new Symbol(workflowExecutionRow.workflowExecutionId.get,
                symbol.key.scope, symbol.key.name, symbol.key.iteration,
                if (symbol.isInput) "INPUT" else "OUTPUT",
                symbol.wdlType.toString, symbol.wdlValue.map(wdlValueToString).orNull)
            }
        }

      }
    }

    val allInsertsOrUpdates = DBIO.sequence(callInsertsOrUpdates.toSeq ++ symbolInsertOrUpdates.toSeq.flatten)
    val allInsertsOrUpdatesFuture = database.run(allInsertsOrUpdates.transactionally)
    Await.result(allInsertsOrUpdatesFuture, Duration.Inf)
  }

  override def updateWorkflow(id: WorkflowId, state: WorkflowState): Unit = {
    val updateFuture: Future[Int] = database.run((for {
      workflowExecutionSelect <- dataAccess.workflowExecutions.filter(_.workflowExecutionUuid === id.toString)
    } yield workflowExecutionSelect.status).update(state.toString).transactionally)
    // TODO: "Assert" only one row updated?
    Await.result(updateFuture, Duration.Inf)
  }

  override def query(workflowId: Option[Seq[WorkflowId]], wdlUris: Option[Seq[String]],
                     states: Option[Seq[WorkflowState]], beforeStart: Option[Date],
                     afterStart: Option[Date], beforeEnd: Option[Date], afterEnd: Option[Date]):
  Seq[QueryWorkflowExecutionResult] = {

    // TODO: Another way to construct a query without an empty filter?
    var query = dataAccess.workflowExecutions.filter(row => true: Rep[Boolean])

    if (workflowId.isDefined) {
      query = query.filter { row =>
        val reps = workflowId.get.map { id =>
          row.workflowExecutionUuid === id.toString
        }
        // Modified from: http://slick.typesafe.com/doc/3.0.0/queries.html#sorting-and-filtering
        reps.reduceLeftOption(_ || _).getOrElse(true: Rep[Boolean])
      }
    }

    if (wdlUris.isDefined) {
      query = query.filter { row =>
        val reps = wdlUris.get.map { id =>
          row.wdlUri === id
        }
        reps.reduceLeftOption(_ || _).getOrElse(true: Rep[Boolean])
      }
    }

    if (states.isDefined) {
      query = query.filter { row =>
        val reps = states.get.map { id =>
          row.status === id.toString
        }
        reps.reduceLeftOption(_ || _).getOrElse(true: Rep[Boolean])
      }
    }

    if (beforeStart.isDefined) {
      query = query.filter(_.startDt < date2timestamp(beforeStart.get))
    }

    if (afterStart.isDefined) {
      query = query.filter(_.startDt > date2timestamp(afterStart.get))
    }

    if (beforeEnd.isDefined) {
      query = query.filter(_.endDt < date2timestamp(beforeEnd.get))
    }

    if (afterEnd.isDefined) {
      query = query.filter(_.endDt > date2timestamp(afterEnd.get))
    }

    // TODO: For now, running queries in separate transactions, and use futures to compose the results.
    // Find more elegant one-to-many query option.
    // In 2.1 we could pass in a transaction session. Not sure how to do this in 3.0.
    // It may be that Slick's composition may just "work" this way.
    // http://slick.typesafe.com/doc/3.0.0/orm-to-slick.html#navigating-the-object-graph
    // http://stackoverflow.com/questions/30567365/slick-3-0-0-how-to-query-one-to-many-many-to-many-relations

    val workflowExecutionRowsFuture = database.run(query.result.transactionally)
    val resultsFuture = workflowExecutionRowsFuture flatMap { workflowExecutionRows =>

      val results = workflowExecutionRows map { workflowExecutionRow =>

        // Query all the symbol store entries
        val symbolsQuery = dataAccess.symbols.filter(
          _.workflowExecutionId === workflowExecutionRow.workflowExecutionId.get)
        val symbolsFuture: Future[Seq[Symbol]] = database.run(symbolsQuery.result.transactionally)

        // Query all the local jobs
        val localJobQuery = for {
          executionRow <- dataAccess.executions.filter(
            _.workflowExecutionId === workflowExecutionRow.workflowExecutionId.get)
          localJobRow <- dataAccess.localJobs.filter(
            _.executionId === executionRow.executionId)
        } yield (executionRow, localJobRow)
        val localJobQueryFuture = database.run(localJobQuery.result.transactionally)

        // Query all the jes jobs
        val jesJobQuery = for {
          executionRow <- dataAccess.executions.filter(
            _.workflowExecutionId === workflowExecutionRow.workflowExecutionId.get)
          jesJobRow <- dataAccess.jesJobs.filter(
            _.executionId === executionRow.executionId)
        } yield (executionRow, jesJobRow)
        val jesJobQueryFuture = database.run(jesJobQuery.result.transactionally)

        val rowTuples = for {
          symbols <- symbolsFuture
          jesJobs <- jesJobQueryFuture
          localJobs <- localJobQueryFuture
        } yield (symbols, jesJobs, localJobs)

        rowTuples map { rowTuple =>

          val symbolStoreEntries = rowTuple._1 map { symbolRow =>
            val wdlType: WdlType = WdlType.fromWdlString(symbolRow.wdlType)
            new SymbolStoreEntry(
              new SymbolStoreKey(
                symbolRow.scope,
                symbolRow.name,
                symbolRow.iteration,
                symbolRow.io == "INPUT"),
              wdlType,
              Option(symbolRow.wdlValue).map(value => wdlType.coerceRawValue(value).get))
          }

          val localJobInfos = rowTuple._3 map { case (executionRow, localJobRow) =>
            new LocalCallInfo(
              executionRow.callFqn,
              ExecutionStatus.withName(executionRow.status),
              localJobRow.pid,
              localJobRow.command,
              localJobRow.rc)
          }

          val jesJobInfos = rowTuple._2 map { case (executionRow, jesJobRow) =>
            new JesCallInfo(
              executionRow.callFqn,
              ExecutionStatus.withName(executionRow.status),
              jesJobRow.jesId,
              jesJobRow.jesStatus)
          }

          new QueryWorkflowExecutionResult(
            UUID.fromString(workflowExecutionRow.workflowExecutionUuid),
            workflowExecutionRow.wdlUri,
            WorkflowState.fromString(workflowExecutionRow.status),
            workflowExecutionRow.startDt,
            workflowExecutionRow.endDt,
            (localJobInfos ++ jesJobInfos).toSet,
            symbolStoreEntries.toSet,
            """{"TODO": "need wdl source"}""", // TODO: Is this in the database schema? If not, add it.
            """{"TODO": "need wdl raw inputs"}""") // TODO: Is this in the database schema? If not, add it.
        }
      }

      // For flatMapping, convert Seq[Future] => Future[Seq]
      Future.sequence(results)
    }

    Await.result(resultsFuture, Duration.Inf)
  }
}
