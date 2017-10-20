package cromwell.core

import wom.graph.GraphNodePort.OutputPort

object CromwellGraphNode {
  // TODO WOM: we could move this into WOM as it's already there for GraphNodes
  implicit class CromwellEnhancedOutputPort(val outputPort: OutputPort) extends AnyVal {
    def unqualifiedName = outputPort.identifier.localName.value
    def fullyQualifiedName = outputPort.identifier.fullyQualifiedName.value
  }
}
