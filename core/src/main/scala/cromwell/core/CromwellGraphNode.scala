package cromwell.core

import wom.callable.Callable.OutputDefinition
import wom.graph.GraphNode
import wom.graph.GraphNodePort.{InputPort, OutputPort}

// TODO WOM: https://github.com/broadinstitute/wdl4s/issues/193
object CromwellGraphNode {
  implicit class CromwellEnhancedGraphNode(val graphNode: GraphNode) extends AnyVal {
    def unqualifiedName = graphNode.name
    def fullyQualifiedName = graphNode.name
  }

  implicit class CromwellEnhancedOutputDefinition(val outputDefinition: OutputDefinition) extends AnyVal {
    def unqualifiedName = outputDefinition.name
    def fullyQualifiedName = outputDefinition.name
  }

  implicit class CromwellEnhancedOutputPort(val outputPort: OutputPort) extends AnyVal {
    def unqualifiedName = outputPort.name
    def fullyQualifiedName = outputPort.graphNode.fullyQualifiedName + "." + outputPort.unqualifiedName
  }

  implicit class CromwellEnhancedInputPort(val inputPort: InputPort) extends AnyVal {
    def unqualifiedName = inputPort.name
    def fullyQualifiedName = inputPort.graphNode.fullyQualifiedName + "." + inputPort.unqualifiedName
  }
}
