package wdl.transforms.draft2.wdlom2wom

import common.validation.ErrorOr.ErrorOr
import wdl.draft2.model.{AstTools, Scope, WdlWorkflow}
import wdl.draft2.parser.WdlParser.Terminal
import wom.callable.{MetaValueElement, WorkflowDefinition}
import wom.transforms.WomWorkflowDefinitionMaker
import wom.transforms.WomGraphMaker.ops._
import wom.SourceFileLocation

object WdlDraft2WomWorkflowDefinitionMaker extends WomWorkflowDefinitionMaker[WdlWorkflow] {
  override def toWomWorkflowDefinition(wdlWorkflow: WdlWorkflow, isASubworkflow: Boolean): ErrorOr[WorkflowDefinition] = {
    // NB: We don't allow "OuterGraphInputNode"s when building this (the Map is empty), so preserveScatterForExternalLookups isn't ever actually used.

    // Figure out the start line of the workflow in the source file
    val t: Terminal = AstTools.findTerminals(wdlWorkflow.ast).head

    def stringifyMetaValues(meta: Map[String, String]): Map[String, MetaValueElement] = {
      meta map {
        case (key, value) =>
          key -> MetaValueElement.MetaValueElementString(value)
      }
    }

    (wdlWorkflow: Scope).toWomGraph(Set.empty, Map.empty, preserveIndexForOuterLookups = true, isASubworkflow: Boolean) map { wg =>
      WorkflowDefinition(
        wdlWorkflow.fullyQualifiedName,
        wg,
        stringifyMetaValues(wdlWorkflow.meta),
        stringifyMetaValues(wdlWorkflow.parameterMeta),
        Some(SourceFileLocation(t.getLine)))
    }
  }
}
