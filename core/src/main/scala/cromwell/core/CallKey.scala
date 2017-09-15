package cromwell.core

import wdl4s.wom.graph.CallNode

trait CallKey extends JobKey {
  def node: CallNode
}
