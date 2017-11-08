package wom.graph

import wom.expression.WomExpression
import wom.graph.GraphNodePort.GraphNodeOutputPort
import wom.types.{WomOptionalType, WomType}

sealed trait GraphInputNode extends GraphNode {
  def womType: WomType
  lazy val singleOutputPort: GraphNodeOutputPort = GraphNodeOutputPort(localName, womType, this)

  override val inputPorts: Set[GraphNodePort.InputPort] = Set.empty
  override val outputPorts: Set[GraphNodePort.OutputPort] = Set(singleOutputPort)
}

sealed trait ExternalGraphInputNode extends GraphInputNode {
  /**
    * The fully qualified name should be the same as the one we expect the key in the input file to have.
    * e.g in WDL:
    * workflow.wdl:
    *   workflow w {
    *     String s # "name" = "s", "fullyQualifiedIdentifier" = "w.s"
    *   }
    * 
    * input.json:
    *   {
    *     "w.s": "hi!"
    *   }
    * 
    * e.g in CWL:
    * workflow.cwl:
    *   class: Workflow
    *   inputs:
    *     s: string # "name" = "s", "fullyQualifiedIdentifier" = "s"
    *   
    * inputs.yml:
    *   s: "hi !"
    * 
    */
  
  override lazy val singleOutputPort: GraphNodeOutputPort = GraphNodeOutputPort(identifier, womType, this)
}

final case class RequiredGraphInputNode(override val identifier: WomIdentifier,
                                        womType: WomType) extends ExternalGraphInputNode

final case class OptionalGraphInputNode(override val identifier: WomIdentifier,
                                        womType: WomOptionalType) extends ExternalGraphInputNode

// If we want to allow defaults to be "complex" expressions with dependencies we may need to make it an InstantiatedExpression here instead
final case class OptionalGraphInputNodeWithDefault(override val identifier: WomIdentifier,
                                                   womType: WomType,
                                                   default: WomExpression) extends ExternalGraphInputNode

object OuterGraphInputNode {
  def apply(forIdentifier: WomIdentifier, linkToOuterGraph: GraphNodePort.OutputPort, preserveScatterIndex: Boolean): OuterGraphInputNode = {
    new OuterGraphInputNode(forIdentifier.copy(fullyQualifiedName = forIdentifier.fullyQualifiedName.prefixWith("^")), linkToOuterGraph, preserveScatterIndex)
  }
}

/**
  * Used to represent an input to any GraphNode's inner graph which is a link to a value somewhere in the outer graph.
  */
class OuterGraphInputNode protected(override val identifier: WomIdentifier, val linkToOuterGraph: GraphNodePort.OutputPort, val preserveScatterIndex: Boolean) extends GraphInputNode {
  override def womType: WomType = linkToOuterGraph.womType
  override lazy val singleOutputPort: GraphNodeOutputPort = GraphNodeOutputPort(identifier, womType, this)
  lazy val linkToOuterGraphNode = linkToOuterGraph.graphNode
}

final case class ScatterVariableNode(override val identifier: WomIdentifier,
                                     scatterExpressionOutputPort: GraphNodePort.OutputPort,
                                     override val womType: WomType) extends OuterGraphInputNode(identifier, scatterExpressionOutputPort, preserveScatterIndex = true)
