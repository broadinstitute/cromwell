package cromwell.engine.db

import java.util.Date

import cromwell.binding.WdlSource
import cromwell.engine.{SymbolStoreEntry, WorkflowId, WorkflowState}

case class QueryWorkflowExecutionResult(workflowId: WorkflowId, wdlUri: String, state: WorkflowState,
                                        startTime: Date, endTime: Option[Date],
                                        calls: Set[CallBackendInfo], symbols: Set[SymbolStoreEntry],
                                        wdlSource: WdlSource, jsonInputs: String)
