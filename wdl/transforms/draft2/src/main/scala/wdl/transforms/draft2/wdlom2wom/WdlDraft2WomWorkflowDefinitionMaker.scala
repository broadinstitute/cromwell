package wdl.transforms.draft2.wdlom2wom

import common.validation.ErrorOr.ErrorOr
import wdl.{Scope, WdlWorkflow}
import wom.callable.WorkflowDefinition
import wom.transforms.WomWorkflowDefinitionMaker
import wom.transforms.WomGraphMaker.ops._

object WdlDraft2WomWorkflowDefinitionMaker extends WomWorkflowDefinitionMaker[WdlWorkflow] {
  override def toWomWorkflowDefinition(wdlWorkflow: WdlWorkflow): ErrorOr[WorkflowDefinition] = {
    // NB: We don't allow "OuterGraphInputNode"s when building this (the Map is empty), so preserveScatterForExternalLookups isn't ever actually used.

    (wdlWorkflow: Scope).toWomGraph(Set.empty, Map.empty, preserveIndexForOuterLookups = true) map { wg =>
      WorkflowDefinition(
        wdlWorkflow.unqualifiedName,
        wg,
        wdlWorkflow.meta,
        wdlWorkflow.parameterMeta,
        List.empty)
    }
  }
}
