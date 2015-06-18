package cromwell.engine.db

import java.util.Date

import cromwell.binding.Call
import cromwell.engine.{SymbolStoreEntry, WorkflowId, WorkflowState}

object DataAccess {
  // Will stamp the start_dt column with DateTime.Now.
  // Will stamp the WorkflowState as Submitted.
  def createWorkflow(id: WorkflowId, wdlUri: String, symbols: Seq[SymbolStoreEntry]): Unit = ???

  def updateWorkflow(id: WorkflowId, state: WorkflowState): Unit = ???

  def updateCall(id: WorkflowId, call: Call, callStatus: Option[CallStatus],
                 callInfo: Option[CallInfo], symbols: Option[Seq[SymbolStoreEntry]]): Unit = ???

  def query(workflowId: Option[Seq[WorkflowId]], wdlUris: Option[Seq[String]], states: Option[Seq[WorkflowState]],
            beforeStart: Option[Date], afterStart: Option[Date], beforeEnd: Option[Date], afterEnd: Option[Date]):
  Seq[QueryWorkflowExecutionResult] = ???
}
