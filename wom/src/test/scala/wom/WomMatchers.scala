package wom

import org.scalactic.Equality
import wom.graph.GraphNode

trait WomMatchers {
  // This will take precedence when comparing graph nodes or collections of graph nodes
  implicit val graphNodeReferenceEquality = new Equality[GraphNode] {
    override def areEqual(left: GraphNode, right: Any): Boolean = right match {
      case node: GraphNode => left eq node
      case _ => false
    }
  }
}

object WomMatchers extends WomMatchers
