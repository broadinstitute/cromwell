package cromwell.engine.db.slick

import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, FullyQualifiedName}
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine.backend.Backend
import cromwell.engine.db.DataAccess.WorkflowInfo
import cromwell.engine.db._
import cromwell.engine.{SymbolStoreKey, SymbolStoreEntry, WorkflowId, WorkflowState}

import scala.concurrent.Future
import scala.language.implicitConversions

// TODO: Need lots of refactoring including: better futures handling, compiled queries, utility functions, etc.
object DataAccessController extends DataAccess {
  /**
   * Creates a row in each of the backend-info specific tables for each call in `calls` corresponding to the backend
   * `backend`.  Or perhaps defer this?
   */
  override def createWorkflow(workflowStateResult: WorkflowInfo, workflowInputs: Traversable[SymbolStoreEntry], calls: Traversable[Call], backend: Backend): Future[Unit] = ???

  override def setStatus(workflowId: WorkflowId, calls: Traversable[Call], callStatus: CallStatus): Future[Unit] = ???

  /** Returns all outputs for this workflowId */
  override def getOutputs(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]] = ???

  /** Get all outputs for the scope of this call. */
  override def getOutputs(workflowId: WorkflowId, call: Call): Future[Traversable[SymbolStoreEntry]] = ???

  /** Should fail if a value is already set.  The keys in the Map are locally qualified names. */
  override def setOutputs(workflowId: WorkflowId, call: Call, callOutputs: Map[String, WdlValue]): Future[Unit] = ???

  override def getExecutionStatuses(workflowId: WorkflowId): Future[Map[FullyQualifiedName, CallStatus]] = ???

  override def getExecutionBackendInfo(workflowId: WorkflowId, call: Call): Future[CallBackendInfo] = ???

  override def updateWorkflowState(workflowId: WorkflowId, workflowState: WorkflowState): Future[Unit] = ???

  override def getWorkflowState(workflowId: WorkflowId): Future[Option[WorkflowState]] = ???

  override def updateExecutionBackendInfo(workflowId: WorkflowId, call: Call, backendInfo: CallBackendInfo): Future[Unit] = ???

  override def getWorkflowsByState(states: Traversable[WorkflowState]): Future[Traversable[WorkflowInfo]] = ???

  /** Get all inputs for the scope of this call. */
  override def getInputs(id: WorkflowId, call: Call): Future[Traversable[SymbolStoreEntry]] = ???
}
