package cromwell.database.slick

import java.sql.{Clob, Timestamp}

import cromwell.database.sql.OldeWorldeSqlDatabase
import cromwell.database.sql.tables.{Execution, ExecutionInfo, RuntimeAttribute, Symbol, WorkflowExecution, WorkflowExecutionAux, WorkflowMetadataSummary}

import scala.concurrent.{ExecutionContext, Future}

trait OldeWorldeSlickDatabase extends OldeWorldeSqlDatabase {
  this: SlickDatabase =>

  import dataAccess.driver.api._

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

  override def runningExecutionsAndExecutionInfos(workflowUuid: String, statuses: Set[String])
                                                 (implicit ec: ExecutionContext):
  Future[Traversable[(Execution, ExecutionInfo)]] = {
    val action = dataAccess.runningExecutionsAndExecutionInfosByWorkflowExecutionUuid(workflowUuid, statuses).result

    runTransaction(action)
  }

}
