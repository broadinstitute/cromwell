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

  protected def setTerminalWithoutClone(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                        statusString: String, endDt: Option[Timestamp], scriptReturnCode: Option[Int],
                                        overallHash: Option[String], dockerHash: Option[String])
                                       (implicit ec: ExecutionContext): Future[Unit]

  protected def setTerminalWithClone(workflowUuid: String, callFqn: String, index: Int, attempt: Int,
                                     statusString: String, endDt: Option[Timestamp], scriptReturnCode: Option[Int],
                                     overallHash: Option[String], dockerHash: Option[String],
                                     workflowUuidClone: String, callFqnClone: String, indexClone: Int,
                                     attemptClone: Int)(implicit ec: ExecutionContext): Future[Unit]

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

  protected def queryWorkflowSummaries(statuses: Set[String], names: Set[String], uuids: Set[String],
                                       startDate: Option[Timestamp], endDate: Option[Timestamp],
                                       page: Option[Int], pageSize: Option[Int])
                                      (implicit ec: ExecutionContext): Future[Traversable[WorkflowMetadataSummary]]

  protected def countWorkflowSummaries(statuses: Set[String], names: Set[String], uuids: Set[String],
                                       startDate: Option[Timestamp], endDate: Option[Timestamp])
                                      (implicit ec: ExecutionContext): Future[Int]

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

  /**
    * Add metadata events to the database transactionally. normalized type structure is as follows:
    * (WorkflowId, MetadataKey, Option[CallFqn, CallIndex, CallAttempt], MetadataValue, MetadataValueType, Timestamp)
    */
  protected def addMetadata(events: Iterable[Metadatum])
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

  /**
    * Retrieves all summarizable metadata satisfying the specified criteria.
    *
    * @param startMetadataId        The minimum ID an entry in `METADATA_JOURNAL` must have to be examined for summary.
    * @param startMetadataTimestamp An optional timestamp.  If specified, a metadatum must have a timestamp greater than or equal to this value.
    * @param buildUpdatedSummary    Takes in the optional existing summary and the metadata, returns the new summary.
    *
    * @return A `Future` with the maximum ID value of the metadata summarized by the invocation of this method.
    */
  protected def refreshMetadataSummaries(startMetadataId: Long, startMetadataTimestamp: Option[Timestamp],
                                         buildUpdatedSummary:
                                         (Option[WorkflowMetadataSummary], Seq[Metadatum]) => WorkflowMetadataSummary)
                                        (implicit ec: ExecutionContext): Future[Long]

  protected def getStatus(workflowUuid: String)
                         (implicit ec: ExecutionContext): Future[Option[String]]
}
