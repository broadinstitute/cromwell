package cromwell.core

import wom.graph.CallNode

trait CallKey extends JobKey {
  def node: CallNode
}
