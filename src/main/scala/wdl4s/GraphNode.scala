package wdl4s

trait GraphNode {
  def upstream: Set[Scope with GraphNode]
  def downstream: Set[Scope with GraphNode]
}
