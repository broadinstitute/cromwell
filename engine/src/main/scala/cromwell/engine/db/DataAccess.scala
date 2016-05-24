package cromwell.engine.db

import cromwell.core.{CallOutputs, WorkflowId}
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine._
import cromwell.engine.backend.{Backend, BackendCallJobDescriptor, _}
import cromwell.engine.db.DataAccess.{ExecutionKeyToJobKey, WorkflowExecutionAndAux}
import cromwell.engine.db.slick._
import cromwell.engine.workflow.{BackendCallKey, ExecutionStoreKey}
import cromwell.webservice.{CallCachingParameters, QueryMetadata, WorkflowQueryParameters, WorkflowQueryResponse}
import wdl4s.values.WdlValue
import wdl4s.{CallInputs, _}

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps

object DataAccess {
  val globalDataAccess: DataAccess = new slick.SlickDataAccess()
  case class ExecutionKeyToJobKey(executionKey: ExecutionDatabaseKey, jobKey: BackendJobKey)
  case class WorkflowExecutionAndAux(execution: WorkflowExecution, aux: WorkflowExecutionAux)
}

trait DataAccess extends AutoCloseable {
  /**
   * Creates a row in each of the backend-info specific tables for each call in `calls` corresponding to the backend
   * `backend`.  Or perhaps defer this?
   */
  def createWorkflow(workflowDescriptor: WorkflowDescriptor,
                     workflowInputs: Traversable[SymbolStoreEntry],
                     calls: Traversable[Scope],
                     backend: Backend)(implicit ec: ExecutionContext): Future[Unit]

  def getWorkflowState(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Option[WorkflowState]]

  def getWorkflowExecutionAndAux(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowExecutionAndAux]

  def getWorkflowExecutionAndAux(workflowExecutionId: Int)(implicit ec: ExecutionContext): Future[WorkflowExecutionAndAux]

  def getWorkflowExecutionAndAuxByState(states: Traversable[WorkflowState])(implicit ec: ExecutionContext): Future[Traversable[WorkflowExecutionAndAux]]

  def getExecutionInfos(workflowId: WorkflowId, call: Call, attempt: Int)(implicit ec: ExecutionContext): Future[Traversable[ExecutionInfo]]

  def updateExecutionInfo(workflowId: WorkflowId, callKey: BackendCallKey, key: String, value: Option[String])
                         (implicit ec: ExecutionContext): Future[Unit]

  def upsertExecutionInfo(workflowId: WorkflowId, callKey: BackendCallKey, keyValues: Map[String, Option[String]])
                         (implicit ec: ExecutionContext): Future[Unit]

  def updateWorkflowState(workflowId: WorkflowId, workflowState: WorkflowState)
                         (implicit ec: ExecutionContext): Future[Unit]

  def getAllSymbolStoreEntries(workflowId: WorkflowId)
                              (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]]

  // TODO needed to support compatibility with current code, this seems like an inefficient way of getting
  // TODO workflow outputs.
  /** Returns all outputs for this workflowId */
  def getWorkflowOutputs(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]]

  def getAllOutputs(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]]

  def getAllInputs(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]]

  /** Get all outputs for the scope of this call. */
  def getOutputs(workflowId: WorkflowId, key: ExecutionDatabaseKey)
                (implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]]

  def upsertRuntimeAttributes(id: WorkflowId, key: ExecutionDatabaseKey, attributes: Map[String, WdlValue])(implicit ec: ExecutionContext): Future[Unit]

  def getAllRuntimeAttributes(id: WorkflowId)(implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, Map[String, String]]]

  /** Get all inputs for the scope of this call. */
  def getInputs(id: WorkflowId, call: Call)(implicit ec: ExecutionContext): Future[Traversable[SymbolStoreEntry]]

  /** Should fail if a value is already set.  The keys in the Map are locally qualified names. */
  def setOutputs(workflowId: WorkflowId, key: ExecutionDatabaseKey, callOutputs: CallOutputs,
                 workflowOutputFqns: Seq[ReportableSymbol])(implicit ec: ExecutionContext): Future[Unit]

  /** Updates the existing input symbols to replace expressions with real values **/
  def updateCallInputs(workflowId: WorkflowId, key: BackendCallKey, callInputs: CallInputs)
                      (implicit ec: ExecutionContext): Future[Int]

  def setExecutionEvents(workflowId: WorkflowId, callFqn: String, shardIndex: Option[Int], attempt: Int,
                         events: Seq[ExecutionEventEntry])(implicit ec: ExecutionContext): Future[Unit]

  /** Gets a mapping from call FQN to an execution event entry list */
  def getAllExecutionEvents(workflowId: WorkflowId)
                           (implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, Seq[ExecutionEventEntry]]]

  /** Add a new failure event for a call into the database. */
  def addCallFailureEvent(workflowId: WorkflowId, executionKey: ExecutionDatabaseKey,
                          failure: FailureEventEntry)(implicit ec: ExecutionContext): Future[Unit]

  /** Add a new failure event for a workflow into the database. */
  def addWorkflowFailureEvent(workflowId: WorkflowId, failure: FailureEventEntry)(implicit ec: ExecutionContext): Future[Unit]

  /** Retrieve all recorded failures for a Workflow */
  def getFailureEvents(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Seq[QualifiedFailureEventEntry]]

  /** Set the status of one or several calls to starting and update the start date. */
  def setStartingStatus(workflowId: WorkflowId, scopeKeys: Traversable[ExecutionDatabaseKey])(implicit ec: ExecutionContext): Future[Unit]

  /** Simply set the status of one of several calls. Status cannot be Starting or a Terminal status. */
  def updateStatus(workflowId: WorkflowId, scopeKeys: Traversable[ExecutionDatabaseKey], status: ExecutionStatus)(implicit ec: ExecutionContext): Future[Unit]

  /** Set the status of a Call to a terminal status, and update associated information (return code, hash, cache). */
  def setTerminalStatus(workflowId: WorkflowId, scopeKeys: ExecutionDatabaseKey, status: ExecutionStatus,
                        scriptReturnCode: Option[Int], hash: Option[ExecutionHash], resultsClonedFrom: Option[BackendCallJobDescriptor])(implicit ec: ExecutionContext): Future[Unit]

  def getExecutionStatuses(workflowId: WorkflowId)
                          (implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, CallStatus]]

  /** Return all execution entries for the FQN, including collector and shards if any */
  def getExecutionStatuses(workflowId: WorkflowId, fqn: FullyQualifiedName)
                          (implicit ec: ExecutionContext): Future[Map[ExecutionDatabaseKey, CallStatus]]

  def getExecutionStatus(workflowId: WorkflowId, key: ExecutionDatabaseKey)
                        (implicit ec: ExecutionContext): Future[Option[CallStatus]]

  def insertCalls(workflowId: WorkflowId, keys: Traversable[ExecutionStoreKey], backend: Backend)
                 (implicit ec: ExecutionContext): Future[Unit]

  def getExecutions(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[Execution]]

  def getExecutionsForRestart(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[Execution]]

  def getExecutionsWithResuableResultsByHash(hash: String)
                                            (implicit ec: ExecutionContext): Future[Traversable[Execution]]

  /** Fetch the workflow having the specified `WorkflowId`. */
  def getWorkflowExecution(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowExecution]

  def getWorkflowExecutionAux(id: WorkflowId)(implicit ec: ExecutionContext): Future[WorkflowExecutionAux]

  def updateWorkflowOptions(workflowId: WorkflowId, workflowOptionsJson: String)(implicit ec: ExecutionContext): Future[Unit]

  def resetTransientExecutions(workflowId: WorkflowId, isTransient: (Execution, Seq[ExecutionInfo]) => Boolean)(implicit ec: ExecutionContext): Future[Unit]

  def findResumableExecutions(workflowId: WorkflowId,
                              isResumable: (Execution, Seq[ExecutionInfo]) => Boolean,
                              jobKeyBuilder: (Execution, Seq[ExecutionInfo]) => BackendJobKey)
                             (implicit ec: ExecutionContext): Future[Traversable[ExecutionKeyToJobKey]]

  def queryWorkflows(queryParameters: WorkflowQueryParameters)
                    (implicit ec: ExecutionContext): Future[(WorkflowQueryResponse, Option[QueryMetadata])]

  def countWorkflows(queryParameters: WorkflowQueryParameters)
                    (implicit ec: ExecutionContext): Future[Int]

  def updateCallCaching(cachingParameters: CallCachingParameters)(implicit ec: ExecutionContext): Future[Int]

  def infosByExecution(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[ExecutionInfosByExecution]]

  def infosByExecution(id: WorkflowId, fqn: FullyQualifiedName)
                      (implicit ec: ExecutionContext): Future[Traversable[ExecutionInfosByExecution]]

  def callCacheDataByExecution(id: WorkflowId)(implicit ec: ExecutionContext): Future[Traversable[ExecutionWithCacheData]]
}
