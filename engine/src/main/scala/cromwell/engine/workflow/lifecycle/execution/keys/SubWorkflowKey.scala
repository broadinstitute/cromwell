package cromwell.engine.workflow.lifecycle.execution.keys

import cromwell.core.CallKey
import cromwell.core.ExecutionIndex._
import wom.graph.WorkflowCallNode

private[execution] case class SubWorkflowKey(node: WorkflowCallNode, index: ExecutionIndex, attempt: Int)
    extends CallKey {
  override val tag = s"SubWorkflow-${node.localName}:${index.fromIndex}:$attempt"
}
