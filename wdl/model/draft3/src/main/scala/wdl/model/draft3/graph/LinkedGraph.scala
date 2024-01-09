package wdl.model.draft3.graph

import wdl.model.draft3.elements.WorkflowGraphElement
import wom.types.WomType

final case class LinkedGraph(elements: Set[WorkflowGraphElement],
                             edges: Set[LinkedGraphEdge],
                             generatedHandles: Set[GeneratedValueHandle],
                             consumedValueLookup: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                             typeAliases: Map[String, WomType]
) {}

final case class LinkedGraphEdge(upstream: WorkflowGraphElement, downstream: WorkflowGraphElement)
