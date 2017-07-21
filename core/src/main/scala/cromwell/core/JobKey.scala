package cromwell.core

import wdl4s.wdl.{Scope, WdlGraphNode}

trait JobKey {
  def scope: Scope with WdlGraphNode
  def index: Option[Int]
  def attempt: Int
  def tag: String

  override def toString = {
    import ExecutionIndex.IndexEnhancedIndex
    s"${scope.fullyQualifiedName}:${index.fromIndex}:$attempt"
  }
  
 def isShard = index.isDefined
}
