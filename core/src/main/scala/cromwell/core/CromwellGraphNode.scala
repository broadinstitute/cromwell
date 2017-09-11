package cromwell.core

import wdl4s.wdl.types.{WdlAnyType, WdlType}
import wdl4s.wom.callable.Callable.OutputDefinition
import wdl4s.wom.graph.GraphNode
import wdl4s.wom.graph.GraphNodePort.OutputPort

object CromwellGraphNode {
  implicit class CromwellEnhancedGraphNode(val graphNode: GraphNode) extends AnyVal {
    def unqualifiedName = graphNode.name.split("#").last
    def fullyQualifiedName = graphNode.name.split("#").last
    // TODO WOM: implement for realz
    def relativeWdlType(node: GraphNode): WdlType = WdlAnyType
  }

  implicit class CromwellEnhancedOutputDefinition(val outputDefinition: OutputDefinition) extends AnyVal {
    def unqualifiedName = outputDefinition.name
    def fullyQualifiedName = outputDefinition.name
    // TODO WOM: implement for realz
    def relativeWdlType(node: GraphNode): WdlType = WdlAnyType
  }

  implicit class CromwellEnhancedOutputPort(val outputPort: OutputPort) extends AnyVal {
    def unqualifiedName = outputPort.name
    def fullyQualifiedName = outputPort.graphNode.fullyQualifiedName + "." + outputPort.unqualifiedName
  }
}
