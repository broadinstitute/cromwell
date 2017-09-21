package cromwell.util

import wom.RuntimeAttributes
import wom.callable.{TaskDefinition, WorkflowDefinition}
import wom.graph.{Graph, TaskCallNode, WorkflowCallNode}

object WomMocks {
  val EmptyTaskDefinition = TaskDefinition("emptyTask", List.empty, RuntimeAttributes(Map.empty),
    Map.empty, Map.empty, List.empty, List.empty)

  val EmptyWorkflowDefinition = mockWorkflowDefinition("emptyWorkflow")

  def mockTaskCall(name: String, definition: TaskDefinition = EmptyTaskDefinition) = {
    TaskCallNode(name, definition, Set.empty, Map.empty)
  }
  
  def mockWorkflowCall(name: String, definition: WorkflowDefinition = EmptyWorkflowDefinition) = {
    WorkflowCallNode(name, definition, Set.empty, Map.empty)
  }

  def mockWorkflowDefinition(name: String) = {
    WorkflowDefinition(name, Graph(Set.empty), Map.empty, Map.empty, List.empty)
  }

  def mockTaskDefinition(name: String) = {
    TaskDefinition(name, List.empty, RuntimeAttributes(Map.empty),
      Map.empty, Map.empty, List.empty, List.empty)
  }
}
