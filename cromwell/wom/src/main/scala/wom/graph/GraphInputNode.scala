package wom.graph

import wom.expression.WomExpression
import wom.graph.GraphNodePort.GraphNodeOutputPort
import wom.graph.expression.ExpressionNode
import wom.types.{WomOptionalType, WomType}

sealed trait GraphInputNode extends GraphNodeWithSingleOutputPort {
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

  /**
    * Key that should be looked for in the input set to satisfy this EGIN
    */
  def nameInInputSet: String
}

final case class RequiredGraphInputNode(override val identifier: WomIdentifier,
                                        womType: WomType,
                                        nameInInputSet: String) extends ExternalGraphInputNode

final case class OptionalGraphInputNode(override val identifier: WomIdentifier,
                                        womType: WomOptionalType,
                                        nameInInputSet: String) extends ExternalGraphInputNode

// If we want to allow defaults to be "complex" expressions with dependencies we may need to make it an InstantiatedExpression here instead
final case class OptionalGraphInputNodeWithDefault(override val identifier: WomIdentifier,
                                                   womType: WomType,
                                                   default: WomExpression,
                                                   nameInInputSet: String) extends ExternalGraphInputNode

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

  lazy val nameToPortMapping: (String, GraphNodeOutputPort) = localName -> singleOutputPort
}

final case class ScatterVariableNode(override val identifier: WomIdentifier,
                                     scatterExpressionNode: ExpressionNode,
                                     override val womType: WomType) extends OuterGraphInputNode(identifier, scatterExpressionNode.singleOutputPort, preserveScatterIndex = true) {
  /*
    * This is the key element of the indexForShard function.
    * Here is an example:
    * Let's say we're scattering over 3 arrays of size 2, 3, and 2 respectively, using a cross product.
    * We will then have 2 * 3 * 2 = 12 shards
    * 
    * The possible combinations of input indices are
    * 
    * Shard #0: (0, 0, 0)
    * Shard #1: (0, 0, 1)
    * Shard #2: (0, 1, 0)
    * Shard #3: (0, 1, 1)
    * Shard #4: (0, 2, 0)
    * Shard #5: (0, 2, 1)
    * Shard #6: (1, 0, 0)
    * Shard #7: (1, 0, 1)
    * Shard #8: (1, 1, 0)
    * Shard #9: (1, 1, 1)
    * Shard #10: (1, 2, 0)
    * Shard #11: (1, 2, 1)
    * 
    * We will also have 3 SVNs, one for each array.
    * Each SVN needs to output one column of the above "matrix", based on the shard number
    * 
    * For SVN 1:          For SVN 2:        For SVN 3:
    *   0 -> 0              0 -> 0            0 -> 0
    *   1 -> 0              1 -> 0            1 -> 1
    *   2 -> 0              2 -> 1            2 -> 0
    *   3 -> 0              3 -> 1            3 -> 1
    *   4 -> 0              4 -> 2            4 -> 0
    *   5 -> 0              5 -> 2            5 -> 1
    *   6 -> 1              6 -> 0            6 -> 0
    *   7 -> 1              7 -> 0            7 -> 1
    *   8 -> 1              8 -> 1            8 -> 0
    *   9 -> 1              9 -> 1            9 -> 1
    *   10 -> 1             10 -> 2           10 -> 0
    *   11 -> 1             11 -> 2           11 -> 1
    *   
    *  shard / 6         (shard / 2) % 3     shard % 2
    *  
    *  This can be generalized to (shard / relativeIndexLength) % arraySize
    *  where 
    *   - shard is the shard index as an Int
    *   - relativeIndexLength is the number of possible combinations for the arrays on the "right" of this variable
    *   - arraySize is the size of the array for this variable
    *   
    * Having a var here is not great as it makes the node mutable, the alternative would be to store this information outside of
    * the node and make it an argument to indexForShard.
   */
  private var _relativeIndexLength: Int = 1
  
  def relativeIndexLength = _relativeIndexLength

  def withRelativeIndexLength(relativeIndexLength: Int) = this._relativeIndexLength = relativeIndexLength

  /**
    * Given a shard number, and the total size of the array, returns at which index in the array the value is located.
    */
  // Maybe this whole function should not be hardcoded here but rather defined by the language.. So far it works for all the use cases we have
  def indexForShard(shardIndex: Int, arraySize: Int) = {
    (shardIndex / _relativeIndexLength) % arraySize
  }
}
