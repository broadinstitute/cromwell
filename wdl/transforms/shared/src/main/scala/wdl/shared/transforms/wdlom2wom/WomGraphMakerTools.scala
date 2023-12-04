package wdl.shared.transforms.wdlom2wom

import wom.graph._
import wom.graph.expression.ExposedExpressionNode

object WomGraphMakerTools {

  // Default outputs should be if we're:
  // - A scatter block
  // - An if block
  // - (NB: top level workflows are already given wildcard outputs in the WDL Workflow building phase)
  def addDefaultOutputs(g: Graph, forWorkflowOutputs: Option[WomIdentifier] = None): Graph = {

    def makeIdentifier(original: WomIdentifier): WomIdentifier = forWorkflowOutputs match {
      case None => original
      case Some(wfIdentifier) => wfIdentifier.combine(original.localName.value)
    }

    Graph(g.nodes.union((g.nodes collect {
      case node: CallNode =>
        node.outputPorts.map { op =>
          val identifier = makeIdentifier(WomIdentifier(op.identifier.localName.value))
          PortBasedGraphOutputNode(identifier, op.womType, op)
        }
      case node: ExposedExpressionNode if forWorkflowOutputs.isEmpty =>
        node.outputPorts.map(op => PortBasedGraphOutputNode(makeIdentifier(WomIdentifier(op.name)), op.womType, op))
      case node: ScatterNode =>
        node.outputMapping.map(op => PortBasedGraphOutputNode(makeIdentifier(op.identifier), op.womType, op))
      case node: ConditionalNode =>
        node.conditionalOutputPorts.map { op =>
          PortBasedGraphOutputNode(makeIdentifier(op.identifier), op.womType, op)
        }
    }).flatten))
  }
}
