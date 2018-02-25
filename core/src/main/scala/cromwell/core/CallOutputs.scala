package cromwell.core

import wom.graph.GraphNodePort.OutputPort
import wom.values.WomValue

object CallOutputs {
  def empty: CallOutputs = CallOutputs(Map.empty)
}
case class CallOutputs(outputs: Map[OutputPort, WomValue])
