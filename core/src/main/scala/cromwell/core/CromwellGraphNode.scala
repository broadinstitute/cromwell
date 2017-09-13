package cromwell.core

import wdl4s.wdl.types.{WdlAnyType, WdlType}
import wdl4s.wom.callable.Callable.OutputDefinition
import wdl4s.wom.graph.GraphNode
import wdl4s.wom.graph.GraphNodePort.{InputPort, OutputPort}

// TODO WOM: https://github.com/broadinstitute/wdl4s/issues/193
object CromwellGraphNode {
  implicit class CromwellEnhancedGraphNode(val graphNode: GraphNode) extends AnyVal {
    def unqualifiedName = graphNode.name
    def fullyQualifiedName = graphNode.name
    def relativeWdlType(node: GraphNode): WdlType = WdlAnyType
  }

  implicit class CromwellEnhancedOutputDefinition(val outputDefinition: OutputDefinition) extends AnyVal {
    def unqualifiedName = outputDefinition.name
    def fullyQualifiedName = outputDefinition.name
    def relativeWdlType(node: GraphNode): WdlType = WdlAnyType
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
