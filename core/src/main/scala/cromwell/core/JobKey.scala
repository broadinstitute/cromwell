package cromwell.core

import wom.graph.GraphNode

trait JobKey {
  def node: GraphNode
  def index: Option[Int]
  def totalIndices: Int = 1
  def attempt: Int
  def tag: String

  override def toString = {
    import ExecutionIndex.IndexEnhancedIndex
    s"${getClass.getSimpleName}_${node.getClass.getSimpleName}_${node.fullyQualifiedName}:${index.fromIndex}:$attempt"
  }
  
 def isShard = index.isDefined
}
