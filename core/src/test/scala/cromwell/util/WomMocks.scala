package cromwell.util

import cats.syntax.validated._
import cromwell.core.CallOutputs
import wom.RuntimeAttributes
import wom.callable.Callable.OutputDefinition
import wom.callable.{CallableTaskDefinition, TaskDefinition, WorkflowDefinition}
import wom.graph.GraphNodePort.{GraphNodeOutputPort, OutputPort}
import wom.graph.{Graph, TaskCallNode, WomIdentifier, WorkflowCallNode}
import wom.types.{WomStringType, WomType}
import wom.values.WomValue

object WomMocks {
  val EmptyTaskDefinition = CallableTaskDefinition("emptyTask", Function.const(List.empty.validNel), RuntimeAttributes(Map.empty),
    Map.empty, Map.empty, List.empty, List.empty, Set.empty, Map.empty)

  val EmptyWorkflowDefinition = mockWorkflowDefinition("emptyWorkflow")

  def mockTaskCall(identifier: WomIdentifier, definition: TaskDefinition = EmptyTaskDefinition) = {
    TaskCallNode(identifier, definition, Set.empty, List.empty)
  }
  
  def mockWorkflowCall(identifier: WomIdentifier, definition: WorkflowDefinition = EmptyWorkflowDefinition) = {
    WorkflowCallNode(identifier, definition, Set.empty, List.empty)
  }

  def mockWorkflowDefinition(name: String) = {
    WorkflowDefinition(name, Graph(Set.empty), Map.empty, Map.empty, List.empty)
  }

  def mockTaskDefinition(name: String) = {
    CallableTaskDefinition(name, Function.const(List.empty.validNel), RuntimeAttributes(Map.empty),
      Map.empty, Map.empty, List.empty, List.empty, Set.empty, Map.empty)
  }
  
  def mockOutputPort(name: String, womType: WomType = WomStringType): OutputPort = {
    GraphNodeOutputPort(WomIdentifier(name, name), womType, null)
  }

  def mockOutputPort(outputDefinition: OutputDefinition): OutputPort = {
    GraphNodeOutputPort(WomIdentifier(outputDefinition.name, outputDefinition.name), outputDefinition.womType, null)
  }
  
  def mockOutputExpectations(outputs: Map[String, WomValue]): CallOutputs = {
    CallOutputs(outputs.map {
      case (key, value) => WomMocks.mockOutputPort(key, value.womType) -> value
    })
  }
}
