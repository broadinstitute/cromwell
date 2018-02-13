package wdl.draft3.transforms.wdlom2wom

import common.Checked
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.ast2wdlom.CheckedAtoB
import wdl.model.draft3.elements.WorkflowDefinitionElement
import wom.callable.WorkflowDefinition
import wom.graph.Graph

object WorkflowDefinitionElementToWomWorkflowDefinition {

  type WorkflowDefinitionElementToWomWorkflowDefinition = CheckedAtoB[WorkflowDefinitionElement, WorkflowDefinition]
  def instance: WorkflowDefinitionElementToWomWorkflowDefinition = CheckedAtoB(convert _)

  def convert(a: WorkflowDefinitionElement): Checked[WorkflowDefinition] = {

    val g: ErrorOr[Graph] = Graph.validateAndConstruct(Set.empty)

    g.map { graph =>
      WorkflowDefinition(a.identifier, graph, Map.empty, Map.empty, List.empty)
    }.toEither

  }
}
