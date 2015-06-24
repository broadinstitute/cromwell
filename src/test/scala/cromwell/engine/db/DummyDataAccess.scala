package cromwell.engine.db

import java.util.Date

import cromwell.binding.{Call, WdlSource}
import cromwell.engine.store.SymbolStoreEntry
import cromwell.engine.{WorkflowId, WorkflowState}

case class DummyDataAccess() extends DataAccess {
  // Will stamp the start_dt column with DateTime.Now.
  override def createWorkflow(id: WorkflowId, wdlUri: String, wdlSource: WdlSource, jsonInputs: String, entries: Seq[SymbolStoreEntry]): Unit = ???

  override def updateCall(id: WorkflowId, call: Call, callStatus: Option[CallStatus], callInfo: Option[CallInfo], symbols: Option[Seq[SymbolStoreEntry]]): Unit = ???

  override def updateWorkflow(id: WorkflowId, state: WorkflowState): Unit = ???

  override def query(workflowId: Option[Seq[WorkflowId]], wdlUris: Option[Seq[String]], states: Option[Seq[WorkflowState]], beforeStart: Option[Date], afterStart: Option[Date], beforeEnd: Option[Date], afterEnd: Option[Date]): Seq[QueryWorkflowExecutionResult] = queryHelper()

  /**
   * None of the arguments to query matter for the current tests, so provide a convenient
   * no-arg method that can be overridden.
   */
  def queryHelper(): Seq[QueryWorkflowExecutionResult] = Seq.empty
}
