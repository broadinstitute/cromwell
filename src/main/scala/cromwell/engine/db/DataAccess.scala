package cromwell.engine.db

import java.util.Date

import cromwell.binding.Call
import cromwell.engine.SymbolStore.SymbolStoreEntry
import cromwell.engine.{WorkflowId, WorkflowState}

trait DataAccess {
  // Will stamp the start_dt column with DateTime.Now.
  // Will stamp the WorkflowState as Submitted.
  def createWorkflow(id: WorkflowId, wdlUri: String, symbols: Seq[SymbolStoreEntry]): Unit

  def updateWorkflow(id: WorkflowId, state: WorkflowState): Unit

  def updateCall(id: WorkflowId, call: Call, callStatus: Option[CallStatus],
                 callInfo: Option[CallInfo], symbols: Option[Seq[SymbolStoreEntry]]): Unit

  def query(workflowId: Option[Seq[WorkflowId]] = None,
            wdlUris: Option[Seq[String]] = None,
            states: Option[Seq[WorkflowState]] = None,
            beforeStart: Option[Date] = None,
            afterStart: Option[Date] = None,
            beforeEnd: Option[Date] = None,
            afterEnd: Option[Date] = None): Seq[QueryWorkflowExecutionResult]
}
