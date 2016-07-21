package cromwell.database

import java.sql.{Clob, SQLTransientException, Timestamp}

import cromwell.database.obj._

import scala.concurrent.{ExecutionContext, Future}
import scalaz.NonEmptyList


trait SqlDatabase extends AutoCloseable with WorkflowStoreSqlDatabase {

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

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  /** Returns all outputs for this workflowId */
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
  def queryWorkflowSummaries(statuses: Set[String], names: Set[String], uuids: Set[String],
                                       startDate: Option[Timestamp], endDate: Option[Timestamp],
                                       page: Option[Int], pageSize: Option[Int])
                                      (implicit ec: ExecutionContext): Future[Traversable[WorkflowMetadataSummary]]

  @deprecated("Olde Worlde Databasee Tablee", "0.21")
  def countWorkflowSummaries(statuses: Set[String], names: Set[String], uuids: Set[String],
                                       startDate: Option[Timestamp], endDate: Option[Timestamp])
                                      (implicit ec: ExecutionContext): Future[Int]

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


  /*
  The following section relates to:
.___  ___.  _______ .___________.    ___       _______       ___   .___________.    ___
|   \/   | |   ____||           |   /   \     |       \     /   \  |           |   /   \
|  \  /  | |  |__   `---|  |----`  /  ^  \    |  .--.  |   /  ^  \ `---|  |----`  /  ^  \
|  |\/|  | |   __|      |  |      /  /_\  \   |  |  |  |  /  /_\  \    |  |      /  /_\  \
|  |  |  | |  |____     |  |     /  _____  \  |  '--'  | /  _____  \   |  |     /  _____  \
|__|  |__| |_______|    |__|    /__/     \__\ |_______/ /__/     \__\  |__|    /__/     \__\
   */

  /**
    * Add metadata events to the database transactionally. normalized type structure is as follows:
    * (WorkflowId, MetadataKey, Option[CallFqn, CallIndex, CallAttempt], MetadataValue, MetadataValueType, Timestamp)
    */
  def addMetadata(events: Iterable[Metadatum])
                           (implicit ec: ExecutionContext): Future[Unit]

  def queryMetadataEvents(workflowUuid: String)
                                   (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  def queryMetadataEvents(workflowUuid: String,
                                    key: String)
                                   (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  def queryMetadataEvents(workflowUuid: String,
                                    callFqn: String,
                                    index: Option[Int],
                                    attempt: Int)
                                   (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  def queryMetadataEvents(workflowUuid: String,
                                    key: String,
                                    callFqn: String,
                                    index: Option[Int],
                                    attempt: Int)
                                   (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  def queryMetadataEventsWithWildcardKeys(workflowUuid: String,
                                                    wildcardKeys: NonEmptyList[String],
                                                    requireEmptyJobKey: Boolean)
                                                   (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  def queryMetadataEventsWithoutWildcardKeys(workflowUuid: String,
                                                       wildcardKeys: NonEmptyList[String],
                                                       requireEmptyJobKey: Boolean)
                                                      (implicit ec: ExecutionContext): Future[Seq[Metadatum]]

  /**
    * Retrieves all summarizable metadata satisfying the specified criteria.
    *
    * @param startMetadataId        The minimum ID an entry in `METADATA_JOURNAL` must have to be examined for summary.
    * @param startMetadataTimestamp An optional timestamp.  If specified, a metadatum must have a timestamp greater than or equal to this value.
    * @param buildUpdatedSummary    Takes in the optional existing summary and the metadata, returns the new summary.
    * @return A `Future` with the maximum ID value of the metadata summarized by the invocation of this method.
    */
  def refreshMetadataSummaries(startMetadataId: Long, startMetadataTimestamp: Option[Timestamp],
                                         buildUpdatedSummary:
                                         (Option[WorkflowMetadataSummary], Seq[Metadatum]) => WorkflowMetadataSummary)
                                        (implicit ec: ExecutionContext): Future[Long]

  def getStatus(workflowUuid: String)
                         (implicit ec: ExecutionContext): Future[Option[String]]

  def workflowExists(possibleWorkflowId: String)
                              (implicit ec: ExecutionContext): Future[Boolean]

}

trait WorkflowStoreSqlDatabase { this: SqlDatabase =>

  def initialize(implicit ec: ExecutionContext): Future[Unit]

  /**
    * Adds the requested WorkflowSourceFiles to the store.
    */
  def add(sources: Iterable[WorkflowStoreEntry])(implicit ec: ExecutionContext): Future[Unit]

  /**
    * Retrieves up to n workflows which have not already been pulled into the engine and sets their pickedUp
    * flag to true
    */
  def fetchRunnableWorkflows(n: Int, state: WorkflowStoreEntryState)(implicit ec: ExecutionContext): Future[List[WorkflowStoreEntry]]

  def remove(id: String)(implicit ec: ExecutionContext): Future[Boolean]
}
