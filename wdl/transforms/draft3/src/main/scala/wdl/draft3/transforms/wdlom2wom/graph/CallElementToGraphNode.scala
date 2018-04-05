package wdl.draft3.transforms.wdlom2wom.graph

import cats.syntax.either._
import cats.instances.list._
import cats.syntax.foldable._
import cats.instances.option._
import cats.instances.vector._
import cats.syntax.validated._
import cats.instances.map._
import cats.syntax.traverse._
import common.validation.ErrorOr.{ErrorOr, _}
import shapeless.Coproduct
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.elements.CallElement
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.callable.Callable.{InputDefinition, InputDefinitionWithDefault, OptionalInputDefinition, RequiredInputDefinition}
import wom.callable.{Callable, CallableTaskDefinition, CommandTaskDefinition}
import wom.graph.CallNode.{CallNodeBuilder, InputDefinitionFold, InputDefinitionPointer}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.expression.{AnonymousExpressionNode, ExpressionNode, PlainAnonymousExpressionNode, TaskCallInputExpressionNode}
import wom.graph._
import wom.types.WomType

object CallElementToGraphNode {
  def convert(a: CallNodeMakerInputs): ErrorOr[Set[GraphNode]] = {
    val callNodeBuilder = new CallNode.CallNodeBuilder()

    // match the call element to a callable
    def callableValidation: ErrorOr[Callable] =
      a.callables.find(_.name == a.node.callableName) match {
        // pass in specific constructor depending on callable type
        case Some(c: Callable) => c.valid
        case None => s"Cannot resolve a callable with name ${a.node.callableName}".invalidNel
      }

    /**
      * Each input definition KV pair becomes an entry in map.
      *
      * i.e.
      * call foo {
      *   input: key = value
      *
      * @return ErrorOr of LocalName(key) mapped to ExpressionNode(value).
      */
    def expressionNodeMappings(callable: Callable): ErrorOr[Map[LocalName, AnonymousExpressionNode]] = {
      a.node.body match {
        case Some(body) =>
          body.inputs.map(input => input.key -> input.value).toMap.traverse {
            case (name, expression) =>
              val identifier = WomIdentifier(name)
              val constructor = callable match {
                case _: CallableTaskDefinition => TaskCallInputExpressionNode.apply _
                case _ => PlainAnonymousExpressionNode.apply _
              }

              AnonymousExpressionNode.fromInputMapping[AnonymousExpressionNode](identifier, WdlomWomExpression(expression, a.linkableValues), a.linkablePorts, constructor) map {
                LocalName(name) -> _
              }
          }
        case None => Map.empty[LocalName, AnonymousExpressionNode].valid
      }
    }

    /**
      * Fold over the input definitions and
      * 1) assign each input definition its InputDefinitionPointer
      * 2) if necessary, create a graph input node and assign its output port to the input definition
      *
      * @return InputDefinitionFold accumulates the input definition mappings, the create graph input nodes, and the expression nodes.
     */
    def foldInputDefinitions(expressionNodes: Map[LocalName, ExpressionNode], callable: Callable): InputDefinitionFold = {
      // Updates the fold with a new graph input node. Happens when an optional or required undefined input without an
      // expression node mapping is found
      def withGraphInputNode(inputDefinition: InputDefinition, graphInputNode: ExternalGraphInputNode) = {
        InputDefinitionFold(
          mappings = List(inputDefinition -> Coproduct[InputDefinitionPointer](graphInputNode.singleOutputPort: OutputPort)),
          callInputPorts = Set(callNodeBuilder.makeInputPort(inputDefinition, graphInputNode.singleOutputPort)),
          newGraphInputNodes = Set(graphInputNode)
        )
      }

      callable.inputs foldMap {
        // If there is an input mapping for this input definition, use that
        case inputDefinition if expressionNodes.contains(inputDefinition.localName) =>
          val expressionNode = expressionNodes(inputDefinition.localName)
          InputDefinitionFold(
            mappings = List(inputDefinition -> expressionNode.inputDefinitionPointer),
            callInputPorts = Set(callNodeBuilder.makeInputPort(inputDefinition, expressionNode.singleOutputPort)),
            newExpressionNodes = Set(expressionNode)
          )

        // No input mapping, use the default expression
        case withDefault@InputDefinitionWithDefault(_, _, expression, _) =>
          InputDefinitionFold(
            mappings = List(withDefault -> Coproduct[InputDefinitionPointer](expression))
          )

        // No input mapping, required and we don't have a default value, create a new RequiredGraphInputNode
        // so that it can be satisfied via workflow inputs
        case required@RequiredInputDefinition(n, womType, _) =>
          val identifier = WomIdentifier(n.value)
          withGraphInputNode(required, RequiredGraphInputNode(identifier, womType, identifier.fullyQualifiedName.value))

        // No input mapping, no default value but optional, create a OptionalGraphInputNode
        // so that it can be satisfied via workflow inputs
        case optional@OptionalInputDefinition(n, womType, _) =>
          val identifier = WomIdentifier(n.value)
          withGraphInputNode(optional, OptionalGraphInputNode(identifier, womType, identifier.fullyQualifiedName.value))
      }
    }

    for {
      callable <- callableValidation
      mappings <- expressionNodeMappings(callable)
      result = callNodeBuilder.build(WomIdentifier(a.node.callableName), callable, foldInputDefinitions(mappings, callable))
    } yield result.nodes
  }
}

case class CallNodeMakerInputs(node: CallElement,
                               linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                               linkablePorts: Map[String, OutputPort],
                               availableTypeAliases: Map[String, WomType],
                               workflowName: String,
                               insideAnotherScatter: Boolean,
                               callables: Set[Callable])
