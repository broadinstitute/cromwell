package wdl.shared.transforms.wdlom2wom

import wom.graph._
import wom.graph.expression.ExposedExpressionNode

object WomGraphMakerTools {

  // Default outputs should be if we're:
  // - A scatter block
  // - An if block
  // - (NB: top level workflows are already given wildcard outputs in the WDL Workflow building phase)
  def addDefaultOutputs(g: Graph, compoundCallIdentifiers: Boolean = true): Graph = Graph(g.nodes.union((g.nodes collect {
    case node: CallNode => node.outputPorts.map(op => {
      val identifier = if (node.compoundOutputIdentifiers) node.identifier.combine(op.name) else WomIdentifier(op.name)
      PortBasedGraphOutputNode(identifier, op.womType, op)
    })
    case node: ExposedExpressionNode => node.outputPorts.map(op => {
      PortBasedGraphOutputNode(WomIdentifier(op.name), op.womType, op)
    })
    case node: ScatterNode => node.outputMapping.map(op => {
      PortBasedGraphOutputNode(op.identifier, op.womType, op)
    })
    case node: ConditionalNode => node.conditionalOutputPorts.map(op => {
      PortBasedGraphOutputNode(op.identifier, op.womType, op)
    })
  }).flatten))
}
