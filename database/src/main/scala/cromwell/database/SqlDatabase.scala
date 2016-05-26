package cromwell.database

import java.sql.{Clob, SQLTransientException, Timestamp}

import cromwell.database.obj._

import scala.concurrent.{ExecutionContext, Future}

trait SqlDatabase extends AutoCloseable {
  protected def isTransient(throwable: Throwable): Boolean = {
    throwable match {
      case e: SQLTransientException => true
      case _ => false
    }
  }

  protected def createWorkflow(workflowExecution: WorkflowExecution,
                               workflowExecutionAuxFromWorkflowExecutionId: Int => WorkflowExecutionAux,
                               workflowInputsFromWorkflowExecutionId: Int => Seq[Symbol],
                               workflowExecutionsFromWorkflowExecutionId: Int => Seq[Execution],
                               workflowExecutionInfosFromExecutionId: Int => Seq[ExecutionInfo])
                              (implicit ec: ExecutionContext): Future[Unit]

  protected def getWorkflowExecution(workflowUuid: String)
                                    (implicit ec: ExecutionContext): Future[WorkflowExecution]

  protected def getWorkflowStateString(workflowUuid: String)(implicit ec: ExecutionContext): Future[Option[String]]

  protected def getWorkflowExecutionAndAuxTuple(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[(WorkflowExecution, WorkflowExecutionAux)]

  protected def getWorkflowExecutionAndAuxTuple(workflowExecutionId: Int)(implicit ec: ExecutionContext):
  Future[(WorkflowExecution, WorkflowExecutionAux)]

  protected def getWorkflowExecutionAndAuxTuples(states: Traversable[String])
                                                (implicit ec: ExecutionContext):
  Future[Traversable[(WorkflowExecution, WorkflowExecutionAux)]]

  protected def getExecutionInfos(workflowUuid: String, callFqn: String, attempt: Int)
                                 (implicit ec: ExecutionContext): Future[Traversable[ExecutionInfo]]

  protected def getExecutionInfoByKey(workflowUuid: String, callFqn: String, attempt: Int, key: String)
                                     (implicit ec: ExecutionContext): Future[Option[Option[String]]]

  protected def updateExecutionInfo(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                    key: String, value: Option[String])(implicit ec: ExecutionContext): Future[Unit]

  protected def upsertExecutionInfo(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                    keyValues: Map[String, Option[String]])(implicit ec: ExecutionContext): Future[Unit]

  protected def updateWorkflowState(workflowUuid: String, workflowState: String, endDate: Option[Timestamp])
                                   (implicit ec: ExecutionContext): Future[Unit]

  protected def getAllSymbols(workflowUuid: String)
                             (implicit ec: ExecutionContext): Future[Traversable[Symbol]]

  /** Returns all outputs for this workflowId */
  protected def getWorkflowOutputSymbols(workflowUuid: String)
                                        (implicit ec: ExecutionContext): Future[Traversable[Symbol]]

  protected def getAllSymbols(workflowUuid: String, ioValue: String)
                             (implicit ec: ExecutionContext): Future[Traversable[Symbol]]

  protected def getAllSymbols(workflowUuid: String, ioValue: String, callFqn: String, index: Int)
                             (implicit ec: ExecutionContext): Future[Traversable[Symbol]]

  protected def setOutputs(workflowUuid: String, symbolsFromWorkflowExecutionId: Int => Seq[Symbol])
                          (implicit ec: ExecutionContext): Future[Unit]

  protected def upsertRuntimeAttributes(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                        attributes: Map[String, String])
                                       (implicit ec: ExecutionContext): Future[Unit]

  protected def getAllRuntimeAttributes(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[Traversable[(String, Int, Int, String, String)]]

  protected def updateCallInputs(workflowUuid: String, callFqn: String, index: Int,
                                 callInputs: Traversable[(String, String, Option[Clob])])
                                (implicit ec: ExecutionContext): Future[Traversable[Int]]

  protected def setExecutionEvents(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                   executionEventsFromExecutionId: Int => Seq[ExecutionEvent])
                                  (implicit ec: ExecutionContext): Future[Unit]

  protected def setExecutionEvents(workflowUuid: String, callFqn: String, attempt: Int,
                                   executionEventsFromExecutionId: Int => Seq[ExecutionEvent])
                                  (implicit ec: ExecutionContext): Future[Unit]

  protected def getAllExecutionEvents(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[Traversable[(String, Int, Int, String, Timestamp, Timestamp)]]

  protected def addCallFailureEvent(workflowUuid: String, callFqn: String, attempt: Int,
                                    failureEventsFromWorkflowExecutionIdAndExecutionId: (Int, Int) => FailureEvent)
                                   (implicit ec: ExecutionContext): Future[Unit]

  protected def addCallFailureEvent(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                    failureEventsFromWorkflowExecutionIdAndExecutionId: (Int, Int) => FailureEvent)
                                   (implicit ec: ExecutionContext): Future[Unit]

  protected def setTerminalWithoutClone(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                        statusString: String, endDt: Option[Timestamp], scriptReturnCode: Option[Int],
                                        overallHash: Option[String], dockerHash: Option[String])
                                       (implicit ec: ExecutionContext): Future[Unit]

  protected def setTerminalWithClone(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                     statusString: String, endDt: Option[Timestamp], scriptReturnCode: Option[Int],
                                     overallHash: Option[String], dockerHash: Option[String],
                                     workflowUuidClone: String, callFqnClone: String, indexClone: Int,
                                     attemptClone: Int)(implicit ec: ExecutionContext): Future[Unit]

  protected def addWorkflowFailureEvent(workflowUuid: String, failureEventsFromWorkflowExecutionId: Int => FailureEvent)
                                       (implicit ec: ExecutionContext): Future[Unit]

  protected def getFailureEvents(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[Traversable[(String, String, Timestamp, Option[String], Option[Int], Option[Int])]]

  protected def setStartingStatus(workflowUuid: String, statusString: String, startDt: Option[Timestamp],
                                  scopeKeys: Traversable[(String, Int, Int)])
                                 (implicit ec: ExecutionContext): Future[Unit]

  protected def updateStatus(workflowUuid: String, statusString: String, scopeKeys: Traversable[(String, Int, Int)])
                            (implicit ec: ExecutionContext): Future[Unit]

  protected def getExecutionStatuses(workflowUuid: String)
                                    (implicit ec: ExecutionContext):
  Future[Traversable[(String, Int, Int, String, Option[Int], Option[String], Option[String])]]

  protected def getExecutionStatuses(workflowUuid: String, callFqn: String)
                                    (implicit ec: ExecutionContext):
  Future[Traversable[(String, Int, Int, String, Option[Int], Option[String], Option[String])]]

  protected def getExecutionStatus(workflowUuid: String, callFqn: String, index: Int, attempt: Int)
                                  (implicit ec: ExecutionContext):
  Future[Option[(String, Option[Int], Option[String], Option[String])]]

  protected def insertCalls(workflowUuid: String,
                            executionsFromWorkflowExecutionId: Int => Seq[Execution],
                            executionInfosFromWorkflowExecutionId: Int => Seq[ExecutionInfo])
                           (implicit ec: ExecutionContext): Future[Unit]

  protected def getExecutions(workflowUuid: String)
                             (implicit ec: ExecutionContext): Future[Traversable[Execution]]

  protected def updateWorkflowOptions(workflowUuid: String, workflowOptionsJson: Clob)
                                     (implicit ec: ExecutionContext): Future[Unit]

  protected def queryWorkflowExecutions(statuses: Set[String], names: Set[String], uuids: Set[String],
                                        startDate: Option[Timestamp], endDate: Option[Timestamp])
                                       (implicit ec: ExecutionContext): Future[Traversable[WorkflowExecution]]

  protected def updateCallCaching(workflowUuid: String, allow: Boolean)(implicit ec: ExecutionContext): Future[Int]

  protected def updateCallCaching(workflowUuid: String, allow: Boolean, callFqn: String, attempt: Int)
                                 (implicit ec: ExecutionContext): Future[Int]

  protected def updateCallCaching(workflowUuid: String, allow: Boolean, callFqn: String, index: Int, attempt: Int)
                                 (implicit ec: ExecutionContext): Future[Int]

  protected def infosByExecution(workflowUuid: String)
                                (implicit ec: ExecutionContext): Future[Traversable[(Execution, ExecutionInfo)]]

  protected def infosByExecution(workflowUuid: String, callFqn: String)
                                (implicit ec: ExecutionContext): Future[Traversable[(Execution, ExecutionInfo)]]

  protected def callCacheDataByExecution(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[Traversable[(Execution, Option[String], Option[String])]]

  protected def getExecutionsWithResuableResultsByHash(hash: String)
                                                      (implicit ec: ExecutionContext): Future[Traversable[Execution]]

  protected def runningExecutionsAndExecutionInfos(workflowUuid: String, statuses: Set[String])
                                                  (implicit ec: ExecutionContext):
  Future[Traversable[(Execution, ExecutionInfo)]]

  protected def addMetadataEvent(workflowUuid: String,
                                 key: String,
                                 value: String,
                                 valueType: String,
                                 timestamp: Timestamp)
                                (implicit ec: ExecutionContext): Future[Unit]

  protected def addMetadataEvent(workflowUuid: String,
                                 key: String,
                                 callFqn: String,
                                 index: Option[Int],
                                 attempt: Int,
                                 value: String,
                                 valueType: String,
                                 timestamp: Timestamp)
                                (implicit ec: ExecutionContext): Future[Unit]

  protected def queryMetadataEvents(workflowUuid: String)
                                   (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  protected def queryMetadataEvents(workflowUuid: String,
                                    key: String)
                                   (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  protected def queryMetadataEvents(workflowUuid: String,
                                    callFqn: String,
                                    index: Option[Int],
                                    attempt: Int)
                                   (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  protected def queryMetadataEvents(workflowUuid: String,
                                    key: String,
                                    callFqn: String,
                                    index: Option[Int],
                                    attempt: Int)
                                   (implicit ec: ExecutionContext): Future[Seq[Metadatum]]
}
