package wdl.draft3.transforms.wdlom2wom

import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.WorkflowDefinitionElement
import wom.callable.WorkflowDefinition
import wom.graph.Graph

object WorkflowDefinitionElementToWomWorkflowDefinition {

  def convert(a: WorkflowDefinitionElement): ErrorOr[WorkflowDefinition] = {

    val g: ErrorOr[Graph] = Graph.validateAndConstruct(Set.empty)

    g.map { graph =>
      WorkflowDefinition(a.identifier, graph, Map.empty, Map.empty, List.empty)
    }

  }
}
