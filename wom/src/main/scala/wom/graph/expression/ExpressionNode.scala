package wom.graph.expression

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.Checked
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}
import common.validation.Validation._
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.GraphNodePort.{ConnectedInputPort, GraphNodeOutputPort, InputPort, OutputPort}
import wom.graph.{GraphNode, GraphNodePort, GraphNodeWithSingleOutputPort, WomIdentifier}
import wom.types.WomType
import wom.values.WomValue

import scala.util.Try

/**
  * Encapsulates a WomExpression with input ports connected to the expression's dependencies.
  */
abstract class ExpressionNode(override val identifier: WomIdentifier,
                              val womExpression: WomExpression,
                              val womType: WomType,
                              val inputMapping: Map[String, InputPort]) extends GraphNodeWithSingleOutputPort with ExpressionNodeLike {
  override val singleOutputPort = GraphNodeOutputPort(identifier, womType, this)
  override val outputPorts: Set[GraphNodePort.OutputPort] = Set(singleOutputPort)
  override val inputPorts = inputMapping.values.toSet
  
  // Again an instance of not so pretty flatMapping with mix of ErrorOrs, Eithers and Tries..
  // TODO: This should return an EitherT or whatever we decide we want to use to package Exceptions + Nel[String]
  /**
    * Evaluates the expression and coerces it.
    */
  def evaluateAndCoerce(inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): Checked[WomValue] = (for {
    evaluated <- womExpression.evaluateValue(inputs, ioFunctionSet)
    coerced <- womType.coerceRawValue(evaluated).toErrorOr
  } yield coerced).leftMap(_.map(e => s"Evaluating ${womExpression.sourceString} failed: $e")).toEither

  override final def evaluate(outputPortLookup: OutputPort => ErrorOr[WomValue], ioFunctionSet: IoFunctionSet): Checked[Map[OutputPort, WomValue]] = {
    import common.validation.ErrorOr._
    for {
      inputs <- inputMapping.traverseValues(inputPort => outputPortLookup(inputPort.upstream)).toEither
      evaluated <- evaluateAndCoerce(inputs, ioFunctionSet)
    } yield Map(singleOutputPort -> evaluated)
  }
}

object ExpressionNode {
  /**
    * Constructs an ExpressionNode or a subclass of an expression node.
    * Note: the WomType is the evaluated type derived from the expression.
    */
  type ExpressionNodeConstructor[E <: ExpressionNode] = (WomIdentifier, WomExpression, WomType, Map[String, InputPort]) => E

  /**
    * Using the passed constructor, attempts to build an expression node from input mappings by linking variable references to other
    * output ports.
    */
  def buildFromConstructor[E <: ExpressionNode](constructor: ExpressionNodeConstructor[E])
                                               (identifier: WomIdentifier, expression: WomExpression, inputMapping: Map[String, OutputPort]): ErrorOr[E] = {
    val graphNodeSetter = new GraphNode.GraphNodeSetter[ExpressionNode]()

    val builtExpressionNode = for {
      combined <- linkWithInputs(graphNodeSetter, expression, inputMapping)
      (evaluatedType, inputPorts) = combined
      expressionNode = constructor(identifier, expression, evaluatedType, inputPorts)
      _ = graphNodeSetter._graphNode = expressionNode
    } yield expressionNode

    def safeSourceString(e: WomExpression) = Try(expression.sourceString).getOrElse("<<expected an expression, none found>>")

    builtExpressionNode.leftMap(_.map(e => s"Cannot build expression for '${identifier.fullyQualifiedName.value} = ${safeSourceString(expression)}': $e"))
  }

  /**
    * Attempts to find an output port for all referenced variables in the expression, and created input ports to connect them together.
    */
  private def linkWithInputs(graphNodeSetter: GraphNode.GraphNodeSetter[ExpressionNode], expression: WomExpression, inputMapping: Map[String, OutputPort]): ErrorOr[(WomType, Map[String, InputPort])] = {
    def linkInput(input: String): ErrorOr[(String, InputPort)] = inputMapping.get(input) match {
      case Some(upstreamPort) => (input, ConnectedInputPort(input, upstreamPort.womType, upstreamPort, graphNodeSetter.get)).validNel
      case None => 
        s"Expression cannot be connected without the input $input (provided: ${inputMapping.toString})".invalidNel
    }

    import common.validation.ErrorOr.ShortCircuitingFlatMap
    for {
      linkedInputList <- expression.inputs.toList traverse linkInput
      linkedInputs = linkedInputList.toMap
      inputTypes = linkedInputs map { case (k, v) => k -> v.womType }
      evaluatedType <- expression.evaluateType(inputTypes)
    } yield (evaluatedType, linkedInputs)
  }
}
