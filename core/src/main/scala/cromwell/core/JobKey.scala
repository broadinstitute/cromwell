package cromwell.core

import cromwell.core.CromwellGraphNode._
import wdl4s.wom.graph.GraphNode

trait JobKey {
  def node: GraphNode
  def index: Option[Int]
  def attempt: Int
  def tag: String

  override def toString = {
    import ExecutionIndex.IndexEnhancedIndex
    s"${node.fullyQualifiedName}:${index.fromIndex}:$attempt"
  }
  
 def isShard = index.isDefined
}
