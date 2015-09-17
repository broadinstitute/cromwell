package cromwell.engine.db

import cromwell.binding._
import cromwell.binding.values.WdlValue
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine.backend.Backend
import cromwell.engine.workflow.{ExecutionStoreKey, OutputKey}
import cromwell.engine.{SymbolStoreEntry, WorkflowId, WorkflowState}

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, Future}

object DataAccess {
  def apply(): DataAccess = new slick.SlickDataAccess()

  /**
   * Creates a DataAccess instance, loans it to the function,
   * and then attempts to shut down the instance.
   *
   * @param f Function to run on the data access.
   * @tparam T Return type of the function.
   * @return Result of calling the function with the dataAccess instance.
   */
  def withDataAccess[T](f: DataAccess => T): T = {
    val dataAccess = DataAccess()
    try {
      f(dataAccess)
    } finally {
      // NOTE: shutdown result thrown away
      Await.ready(dataAccess.shutdown(), Duration.Inf)
    }
  }
}

trait DataAccess {
  /**
   * Creates a row in each of the backend-info specific tables for each call in `calls` corresponding to the backend
   * `backend`.  Or perhaps defer this?
   */
  def createWorkflow(workflowDescriptor: WorkflowDescriptor,
                     workflowInputs: Traversable[SymbolStoreEntry],
                     calls: Traversable[Scope],
                     backend: Backend): Future[Unit]

  def getWorkflowState(workflowId: WorkflowId): Future[Option[WorkflowState]]

  def getWorkflow(workflowId: WorkflowId): Future[WorkflowDescriptor]

  def getWorkflowsByState(states: Traversable[WorkflowState]): Future[Traversable[WorkflowDescriptor]]

  def getExecutionBackendInfo(workflowId: WorkflowId, call: Call): Future[CallBackendInfo]

  def updateExecutionBackendInfo(workflowId: WorkflowId, call: Call, backendInfo: CallBackendInfo): Future[Unit]

  def updateWorkflowState(workflowId: WorkflowId, workflowState: WorkflowState): Future[Unit]

  def getFullyQualifiedName(workflowId: WorkflowId, fqn: FullyQualifiedName): Future[Traversable[SymbolStoreEntry]]

  def getAll(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]]

  // TODO needed to support compatibility with current code, this seems like an inefficient way of getting
  // TODO workflow outputs.
  /** Returns all outputs for this workflowId */
  def getOutputs(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]]

  /** Get all outputs for the scope of this call. */
  def getOutputs(workflowId: WorkflowId, key: ExecutionDatabaseKey): Future[Traversable[SymbolStoreEntry]]

  /** Get all inputs for the scope of this call. */
  def getInputs(id: WorkflowId, call: Call): Future[Traversable[SymbolStoreEntry]]

  /** Should fail if a value is already set.  The keys in the Map are locally qualified names. */
  def setOutputs(workflowId: WorkflowId, key: OutputKey, callOutputs: Map[String, WdlValue]): Future[Unit]

  def setStatus(workflowId: WorkflowId, keys: Traversable[ExecutionDatabaseKey], executionStatus: ExecutionStatus): Future[Unit] = {
    setStatus(workflowId, keys, CallStatus(executionStatus, None))
  }

  def setStatus(workflowId: WorkflowId, keys: Traversable[ExecutionDatabaseKey], callStatus: CallStatus): Future[Unit]

  def getExecutionStatuses(workflowId: WorkflowId): Future[Map[ExecutionDatabaseKey, CallStatus]]

  /** Return all execution entries for the FQN, including collector and shards if any */
  def getExecutionStatuses(workflowId: WorkflowId, fqn: FullyQualifiedName): Future[Map[ExecutionDatabaseKey, CallStatus]]

  def getExecutionStatus(workflowId: WorkflowId, key: ExecutionDatabaseKey): Future[Option[CallStatus]]

  def insertCalls(workflowId: WorkflowId, keys: Traversable[ExecutionStoreKey], backend: Backend): Future[Unit]

  /** Shutdown. NOTE: Should (internally or explicitly) use AsyncExecutor.shutdownExecutor. */
  def shutdown(): Future[Unit]
}
