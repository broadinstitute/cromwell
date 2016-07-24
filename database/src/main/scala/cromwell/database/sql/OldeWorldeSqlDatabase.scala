package cromwell.database.sql

import java.sql.{Clob, Timestamp}

import cromwell.database.sql.tables.{Execution, ExecutionInfo, Symbol, WorkflowExecution, WorkflowExecutionAux, WorkflowMetadataSummary}

import scala.concurrent.{ExecutionContext, Future}

trait OldeWorldeSqlDatabase {
  this: SqlDatabase =>

  /*
  The following section relates to:
  ______    __       _______   _______
 /  __  \  |  |     |       \ |   ____|
|  |  |  | |  |     |  .--.  ||  |__
|  |  |  | |  |     |  |  |  ||   __|
|  `--'  | |  `----.|  '--'  ||  |____
 \______/  |_______||_______/ |_______|

____    __    ____  ______   .______       __       _______   _______
\   \  /  \  /   / /  __  \  |   _  \     |  |     |       \ |   ____|
 \   \/    \/   / |  |  |  | |  |_)  |    |  |     |  .--.  ||  |__
  \            /  |  |  |  | |      /     |  |     |  |  |  ||   __|
   \    /\    /   |  `--'  | |  |\  \----.|  `----.|  '--'  ||  |____
    \__/  \__/     \______/  | _| `._____||_______||_______/ |_______|
 */

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def createWorkflow(workflowExecution: WorkflowExecution,
                     workflowExecutionAuxFromWorkflowExecutionId: Int => WorkflowExecutionAux,
                     workflowInputsFromWorkflowExecutionId: Int => Seq[Symbol],
                     workflowExecutionsFromWorkflowExecutionId: Int => Seq[Execution],
                     workflowExecutionInfosFromExecutionId: Int => Seq[ExecutionInfo])
                    (implicit ec: ExecutionContext): Future[Unit]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getWorkflowExecution(workflowUuid: String)
                          (implicit ec: ExecutionContext): Future[WorkflowExecution]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getWorkflowStateString(workflowUuid: String)(implicit ec: ExecutionContext): Future[Option[String]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getWorkflowExecutionAndAuxTuple(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[(WorkflowExecution, WorkflowExecutionAux)]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getWorkflowExecutionAndAuxTuple(workflowExecutionId: Int)(implicit ec: ExecutionContext):
  Future[(WorkflowExecution, WorkflowExecutionAux)]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getWorkflowExecutionAndAuxTuples(states: Traversable[String])
                                      (implicit ec: ExecutionContext):
  Future[Traversable[(WorkflowExecution, WorkflowExecutionAux)]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getExecutionInfos(workflowUuid: String, callFqn: String, attempt: Int)
                       (implicit ec: ExecutionContext): Future[Traversable[ExecutionInfo]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getExecutionInfoByKey(workflowUuid: String, callFqn: String, attempt: Int, key: String)
                           (implicit ec: ExecutionContext): Future[Option[Option[String]]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def updateExecutionInfo(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                          key: String, value: Option[String])(implicit ec: ExecutionContext): Future[Unit]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def upsertExecutionInfo(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                          keyValues: Map[String, Option[String]])(implicit ec: ExecutionContext): Future[Unit]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def updateWorkflowState(workflowUuid: String, workflowState: String, endDate: Option[Timestamp])
                         (implicit ec: ExecutionContext): Future[Unit]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getAllSymbols(workflowUuid: String)
                   (implicit ec: ExecutionContext): Future[Traversable[Symbol]]

  /** Returns all outputs for this workflowId */
  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getWorkflowOutputSymbols(workflowUuid: String)
                              (implicit ec: ExecutionContext): Future[Traversable[Symbol]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getAllSymbols(workflowUuid: String, ioValue: String)
                   (implicit ec: ExecutionContext): Future[Traversable[Symbol]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getAllSymbols(workflowUuid: String, ioValue: String, callFqn: String, index: Int)
                   (implicit ec: ExecutionContext): Future[Traversable[Symbol]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def setOutputs(workflowUuid: String, symbolsFromWorkflowExecutionId: Int => Seq[Symbol])
                (implicit ec: ExecutionContext): Future[Unit]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def upsertRuntimeAttributes(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                              attributes: Map[String, String])
                             (implicit ec: ExecutionContext): Future[Unit]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getAllRuntimeAttributes(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[Traversable[(String, Int, Int, String, String)]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def updateCallInputs(workflowUuid: String, callFqn: String, index: Int,
                       callInputs: Traversable[(String, String, Option[Clob])])
                      (implicit ec: ExecutionContext): Future[Traversable[Int]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def setTerminalWithoutClone(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                              statusString: String, endDt: Option[Timestamp], scriptReturnCode: Option[Int],
                              overallHash: Option[String], dockerHash: Option[String])
                             (implicit ec: ExecutionContext): Future[Unit]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def setTerminalWithClone(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                           statusString: String, endDt: Option[Timestamp], scriptReturnCode: Option[Int],
                           overallHash: Option[String], dockerHash: Option[String],
                           workflowUuidClone: String, callFqnClone: String, indexClone: Int,
                           attemptClone: Int)(implicit ec: ExecutionContext): Future[Unit]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def setStartingStatus(workflowUuid: String, statusString: String, startDt: Option[Timestamp],
                        scopeKeys: Traversable[(String, Int, Int)])
                       (implicit ec: ExecutionContext): Future[Unit]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def updateStatus(workflowUuid: String, statusString: String, scopeKeys: Traversable[(String, Int, Int)])
                  (implicit ec: ExecutionContext): Future[Unit]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getExecutionStatuses(workflowUuid: String)
                          (implicit ec: ExecutionContext):
  Future[Traversable[(String, Int, Int, String, Option[Int], Option[String], Option[String])]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getExecutionStatuses(workflowUuid: String, callFqn: String)
                          (implicit ec: ExecutionContext):
  Future[Traversable[(String, Int, Int, String, Option[Int], Option[String], Option[String])]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getExecutionStatus(workflowUuid: String, callFqn: String, index: Int, attempt: Int)
                        (implicit ec: ExecutionContext):
  Future[Option[(String, Option[Int], Option[String], Option[String])]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def insertCalls(workflowUuid: String,
                  executionsFromWorkflowExecutionId: Int => Seq[Execution],
                  executionInfosFromWorkflowExecutionId: Int => Seq[ExecutionInfo])
                 (implicit ec: ExecutionContext): Future[Unit]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getExecutions(workflowUuid: String)
                   (implicit ec: ExecutionContext): Future[Traversable[Execution]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def updateWorkflowOptions(workflowUuid: String, workflowOptionsJson: Clob)
                           (implicit ec: ExecutionContext): Future[Unit]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def updateCallCaching(workflowUuid: String, allow: Boolean)(implicit ec: ExecutionContext): Future[Int]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def updateCallCaching(workflowUuid: String, allow: Boolean, callFqn: String, attempt: Int)
                       (implicit ec: ExecutionContext): Future[Int]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def updateCallCaching(workflowUuid: String, allow: Boolean, callFqn: String, index: Int, attempt: Int)
                       (implicit ec: ExecutionContext): Future[Int]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def infosByExecution(workflowUuid: String)
                      (implicit ec: ExecutionContext): Future[Traversable[(Execution, ExecutionInfo)]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def infosByExecution(workflowUuid: String, callFqn: String)
                      (implicit ec: ExecutionContext): Future[Traversable[(Execution, ExecutionInfo)]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def callCacheDataByExecution(workflowUuid: String)(implicit ec: ExecutionContext):
  Future[Traversable[(Execution, Option[String], Option[String])]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def getExecutionsWithResuableResultsByHash(hash: String)
                                            (implicit ec: ExecutionContext): Future[Traversable[Execution]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def runningExecutionsAndExecutionInfos(workflowUuid: String, statuses: Set[String])
                                        (implicit ec: ExecutionContext): Future[Traversable[(Execution, ExecutionInfo)]]
}
