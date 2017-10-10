package cromwell.engine.workflow.lifecycle.execution.keys

import cromwell.core.ExecutionIndex._
import cromwell.core.JobKey
import wom.graph.ExpressionNode
import wom.graph.GraphNodePort.OutputPort

/**
  * Key for expression nodes.
  */
private [execution] case class ExpressionKey(node: ExpressionNode, index: ExecutionIndex) extends JobKey {
  override val attempt = 1
  override lazy val tag = s"Expression-${node.localName}:${index.fromIndex}:$attempt"

  val womType = node.womType
  val singleOutputPort: OutputPort = node.singleExpressionOutputPort

  lazy val upstreamPorts: Map[String, OutputPort] = node.inputMapping map {
    case (key, input) => key -> input.upstream
  }
}
