package wdl4s.wom.graph

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import cats.data.Validated.Valid
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wdl.types.WdlType
import wdl4s.wom.expression.{NamedExpressionInput, WomExpression}
import wdl4s.wom.graph.GraphNodePort.{ConnectedInputPort, InputPort, OutputPort}

sealed trait GraphOutputNode extends GraphNode {
  def name: String
  def womType: WdlType

  final override def outputPorts: Set[GraphNodePort.OutputPort] = Set.empty
}

/**
  * Exposes an existing output port as a graph output.
  */
final case class PortBasedGraphOutputNode(name: String, womType: WdlType, source: OutputPort) extends GraphOutputNode {
  override val inputPorts: Set[GraphNodePort.InputPort] = Set(ConnectedInputPort(name, womType, source, _ => this))
}

/**
  * A graph output which is produced by evaluating an expression.
  *
  * NB: Construct this via ExpressionBasedGraphOutputNode.linkWithInputs(...)
  */
final case class ExpressionBasedGraphOutputNode private(name: String, expression: WomExpression, inputMapping: Map[NamedExpressionInput, InputPort]) extends GraphOutputNode {
  override val womType = expression.womType
  override val inputPorts = inputMapping.values.toSet
}

object ExpressionBasedGraphOutputNode {
  def linkWithInputs(name: String, expression: WomExpression, inputMapping: Map[String, OutputPort]): ErrorOr[ExpressionBasedGraphOutputNode] = {
    val graphNodeSetter = new GraphNode.GraphNodeSetter()

    def linkInput(input: NamedExpressionInput): ErrorOr[(NamedExpressionInput, InputPort)] =
      if (inputMapping.contains(input.name)) {
        Valid((input, ConnectedInputPort(input.name, input.womType, inputMapping(input.name), graphNodeSetter.get)))
      } else {
        s"Expression $name cannot be connected without the input $input.".invalidNel
      }

    val linkedInputs: ErrorOr[Map[NamedExpressionInput, InputPort]] = (expression.inputs.toList traverse linkInput).map(_.toMap)
    val graphOutputNode: ErrorOr[ExpressionBasedGraphOutputNode] = linkedInputs map {li => ExpressionBasedGraphOutputNode(name, expression, li)}

    graphOutputNode.foreach(graphNodeSetter._graphNode = _)
    graphOutputNode
  }
}

