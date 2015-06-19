package cromwell.engine.db

import java.util.Date

import cromwell.binding.WdlSource
import cromwell.engine.SymbolStore.SymbolStoreEntry
import cromwell.engine.{WorkflowId, WorkflowState}

case class QueryWorkflowExecutionResult(workflowId: WorkflowId, wdlUri: String, state: WorkflowState,
                                        startTime: Date, endTime: Option[Date],
                                        calls: Set[CallInfo], symbols: Set[SymbolStoreEntry],
                                        wdlSource: WdlSource, wdlRawInputs: String)
