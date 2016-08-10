package cromwell.core

import wdl4s.Scope

trait JobKey {
  def scope: Scope
  def index: Option[Int]
  def attempt: Int
  def tag: String

  override def toString = {
    import ExecutionIndex.IndexEnhancedIndex
    s"${scope.fullyQualifiedName}:${index.fromIndex}:$attempt"
  }
}
