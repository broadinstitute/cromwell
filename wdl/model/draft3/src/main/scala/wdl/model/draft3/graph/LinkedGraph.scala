package wdl.model.draft3.graph

import wdl.model.draft3.elements.WorkflowGraphElement

// TODO: Most of these fields seem like they'll be useful but haven't got any real use case yet. Double check at some point
final case class LinkedGraph(elements: Set[WorkflowGraphElement],
                             generatedValuesByGraphElement: Map[WorkflowGraphElement, Set[GeneratedValueHandle]],
                             edges: Set[LinkedGraphEdge],
                             consumedValuesByGraphElement: Map[WorkflowGraphElement, Set[UnlinkedConsumedValueHook]],
                             graphElementByGeneratedValueHandle: Map[GeneratedValueHandle, WorkflowGraphElement],
                             consumedValueLookup: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]) {
}

final case class LinkedGraphEdge(upstream: WorkflowGraphElement, downstream: WorkflowGraphElement)
