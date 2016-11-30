package wdl4s

trait GraphNode extends Scope {
  def upstream: Set[GraphNode]
  def downstream: Set[GraphNode]
}
