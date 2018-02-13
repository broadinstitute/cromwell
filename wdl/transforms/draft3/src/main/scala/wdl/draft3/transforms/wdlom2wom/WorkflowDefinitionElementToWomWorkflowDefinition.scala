package wdl.draft3.transforms.wdlom2wom

import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.ast2wdlom.FromAtoB
import wdl.model.draft3.elements.WorkflowDefinitionElement
import wom.callable.WorkflowDefinition
import wom.graph.Graph

case class WorkflowDefinitionElementToWomWorkflowDefinition(inputs: Option[String]) extends FromAtoB[WorkflowDefinitionElement, WorkflowDefinition] {
  override def convert(a: WorkflowDefinitionElement): ErrorOr[WorkflowDefinition] = {

    val g: ErrorOr[Graph] = Graph.validateAndConstruct(Set.empty)

    g map { graph =>
      WorkflowDefinition(a.identifier, graph, Map.empty, Map.empty, List.empty)
    }

  }
}
