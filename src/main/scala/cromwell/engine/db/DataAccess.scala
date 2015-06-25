package cromwell.engine.db

import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, FullyQualifiedName, WdlJson, WdlSource}
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine.backend.Backend
import cromwell.engine.{SymbolStoreKey, SymbolStoreEntry, WorkflowId, WorkflowState}

import scala.concurrent.Future

object DataAccess {
  // TODO PLEASE RENAME ME
  case class WorkflowInfo(workflowId: WorkflowId, wdlSource: WdlSource, wdlJson: WdlJson)
}

trait DataAccess {

  import DataAccess._
  /**
   * Creates a row in each of the backend-info specific tables for each call in `calls` corresponding to the backend
   * `backend`.  Or perhaps defer this?
   */
  def createWorkflow(workflowInfo: WorkflowInfo,
                     workflowInputs: Traversable[SymbolStoreEntry],
                     calls: Traversable[Call],
                     backend: Backend): Future[Unit]

  def getWorkflowState(workflowId: WorkflowId): Future[Option[WorkflowState]]

  def getWorkflowsByState(states: Traversable[WorkflowState]): Future[Traversable[WorkflowInfo]]

  def getExecutionBackendInfo(workflowId: WorkflowId, call: Call): Future[CallBackendInfo]

  def updateExecutionBackendInfo(workflowId: WorkflowId, call: Call, backendInfo: CallBackendInfo): Future[Unit]

  def updateWorkflowState(workflowId: WorkflowId, workflowState: WorkflowState): Future[Unit]

  // TODO needed to support compatibility with current code, this seems like an inefficient way of getting
  // TODO workflow outputs.
  /** Returns all outputs for this workflowId */
  def getOutputs(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]]

  /** Get all outputs for the scope of this call. */
  def getOutputs(workflowId: WorkflowId, call: Call): Future[Traversable[SymbolStoreEntry]]

  /** Get all inputs for the scope of this call. */
  def getInputs(id: WorkflowId, call: Call): Future[Traversable[SymbolStoreEntry]]

  /** Should fail if a value is already set.  The keys in the Map are locally qualified names. */
  def setOutputs(workflowId: WorkflowId, call: Call, callOutputs: Map[String, WdlValue]): Future[Unit]

  def setStatus(workflowId: WorkflowId, calls: Traversable[Call], callStatus: CallStatus): Future[Unit]

  def getExecutionStatuses(workflowId: WorkflowId): Future[Map[FullyQualifiedName, CallStatus]]

}