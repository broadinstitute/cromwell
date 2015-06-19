package cromwell.engine.db

import java.util.Date

import cromwell.binding.Call
import cromwell.engine.SymbolStore.SymbolStoreEntry
import cromwell.engine.{WorkflowId, WorkflowState}

case object DummyDataAccess extends DataAccess {
  // Will stamp the start_dt column with DateTime.Now.
  override def createWorkflow(id: WorkflowId, wdlUri: String, symbols: Seq[SymbolStoreEntry]): Unit = ???

  override def updateCall(id: WorkflowId, call: Call, callStatus: Option[CallStatus], callInfo: Option[CallInfo], symbols: Option[Seq[SymbolStoreEntry]]): Unit = ???

  override def updateWorkflow(id: WorkflowId, state: WorkflowState): Unit = ???

  override def query(workflowId: Option[Seq[WorkflowId]], wdlUris: Option[Seq[String]], states: Option[Seq[WorkflowState]], beforeStart: Option[Date], afterStart: Option[Date], beforeEnd: Option[Date], afterEnd: Option[Date]): Seq[QueryWorkflowExecutionResult] = Seq.empty
}
