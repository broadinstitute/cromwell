package cromwell.util

import wom.RuntimeAttributes
import wom.callable.{CallableTaskDefinition, TaskDefinition, WorkflowDefinition}
import wom.graph.{Graph, TaskCallNode, WomIdentifier, WorkflowCallNode}

object WomMocks {
  val EmptyTaskDefinition = CallableTaskDefinition("emptyTask", List.empty, RuntimeAttributes(Map.empty),
    Map.empty, Map.empty, List.empty, List.empty)

  val EmptyWorkflowDefinition = mockWorkflowDefinition("emptyWorkflow")

  def mockTaskCall(identifier: WomIdentifier, definition: TaskDefinition = EmptyTaskDefinition) = {
    TaskCallNode(identifier, definition, Set.empty, Map.empty)
  }
  
  def mockWorkflowCall(identifier: WomIdentifier, definition: WorkflowDefinition = EmptyWorkflowDefinition) = {
    WorkflowCallNode(identifier, definition, Set.empty, Map.empty)
  }

  def mockWorkflowDefinition(name: String) = {
    WorkflowDefinition(name, Graph(Set.empty), Map.empty, Map.empty, List.empty)
  }

  def mockTaskDefinition(name: String) = {
    CallableTaskDefinition(name, List.empty, RuntimeAttributes(Map.empty),
      Map.empty, Map.empty, List.empty, List.empty)
  }
}
