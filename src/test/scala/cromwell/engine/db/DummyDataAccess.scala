package cromwell.engine.db

import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, FullyQualifiedName}
import cromwell.engine.ExecutionStatus.ExecutionStatus
import cromwell.engine._
import cromwell.engine.backend.Backend
import cromwell.engine.db.DataAccess.WorkflowInfo

import scala.collection.concurrent.TrieMap
import scala.concurrent.Future

case class DummyDataAccess() extends DataAccess {

  private val workflowStates: TrieMap[WorkflowId, WorkflowState] = TrieMap.empty

  val executionStatuses: TrieMap[WorkflowId, TrieMap[String, ExecutionStatus]] = TrieMap.empty

  val symbolStore: TrieMap[WorkflowId, TrieMap[SymbolStoreKey, SymbolStoreEntry]] = TrieMap.empty

  /**
   * Creates a row in each of the backend-info specific tables for each call in `calls` corresponding to the backend
   * `backend`.  Or perhaps defer this?
   * FIXME does not do what the comment above says, fix behavior or comment.
   */
  override def createWorkflow(workflowInfo: WorkflowInfo, symbols: Traversable[SymbolStoreEntry],
                              calls: Traversable[Call], backend: Backend): Future[Unit] = {
    Future.successful {
      val id = workflowInfo.workflowId
      executionStatuses += (id -> TrieMap.empty)
      symbolStore += (id -> TrieMap.empty)
      setStatus(id, calls map { _.fullyQualifiedName }, ExecutionStatus.NotStarted)
      symbols foreach { symbol =>
        symbolStore(id)(symbol.key) = symbol
      }
    }
  }

  override def getWorkflowsByState(states: Traversable[WorkflowState]): Future[Traversable[WorkflowInfo]] = {
    val statesSet = states.toSet
    Future.successful(workflowStates.collect { case (id, state) if statesSet.contains(state) => WorkflowInfo(id, "", "")})
  }

  override def setStatus(workflowId: WorkflowId, callFqns: Traversable[FullyQualifiedName], callStatus: CallStatus): Future[Unit] = {
    Future.successful {
      callFqns foreach { callFqn =>
        executionStatuses(workflowId)(callFqn) = callStatus
      }
    }
  }

  private def getSymbols(workflowId: WorkflowId, inputs: Boolean, scope: Option[String]): Traversable[SymbolStoreEntry] = {
    def passesFilter(key: SymbolStoreKey): Boolean = key.input == inputs && (scope.isEmpty || scope.get == key.scope)
    symbolStore(workflowId).collect { case (key, value) if passesFilter(key) => value }
  }

  /** Returns all outputs for this workflowId */
  override def getOutputs(workflowId: WorkflowId): Future[Traversable[SymbolStoreEntry]] = {
    Future.successful { getSymbols(workflowId, inputs = false, scope = None) }
  }

  /** Get all outputs for the scope of this call. */
  override def getOutputs(workflowId: WorkflowId, call: Call): Future[Traversable[SymbolStoreEntry]] = {
    Future.successful { getSymbols(workflowId, inputs = false, scope = Some(call.fullyQualifiedName)) }
  }

  /** Get all inputs for the scope of this call.  TODO refactor with above. */
  override def getInputs(workflowId: WorkflowId, call: Call): Future[Traversable[SymbolStoreEntry]] = {
    println("looking up inputs for " + call.name)
    Future.successful { getSymbols(workflowId, inputs = true, scope = Some(call.fullyQualifiedName)) }
  }

  /** The keys in the Map are locally qualified names. */
  override def setOutputs(workflowId: WorkflowId, call: Call, callOutputs: Map[String, WdlValue]): Future[Unit] = {
    Future.successful {
      callOutputs foreach { case (name, wdlValue) =>
        val entry = SymbolStoreEntry(call.fullyQualifiedName + "." + name, wdlValue, input = false)
        symbolStore(workflowId)(entry.key) = entry
      }
    }
  }

  override def getExecutionStatuses(workflowId: WorkflowId): Future[Map[FullyQualifiedName, CallStatus]] = {
    Future.successful(executionStatuses(workflowId).toMap)
  }

  override def getExecutionBackendInfo(workflowId: WorkflowId, call: Call): Future[CallBackendInfo] = ???

  override def updateWorkflowState(workflowId: WorkflowId, workflowState: WorkflowState): Future[Unit] = {
    Future.successful(workflowStates += (workflowId -> workflowState))
  }

  override def getWorkflowState(workflowId: WorkflowId): Future[Option[WorkflowState]] = {
    Future.successful(workflowStates.get(workflowId))
  }

  override def updateExecutionBackendInfo(workflowId: WorkflowId, call: Call, backendInfo: CallBackendInfo): Future[Unit] = ???

}
