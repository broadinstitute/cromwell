package wdl.draft3.transforms.wdlom2wom.graph

import cats.instances.list._
import cats.syntax.apply._
import cats.syntax.foldable._
import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, _}
import shapeless.Coproduct
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.elements.CallElement
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.callable.Callable._
import wom.callable.{Callable, CallableTaskDefinition, TaskDefinition, WorkflowDefinition}
import wom.graph.CallNode.{CallNodeAndNewNodes, InputDefinitionFold, InputDefinitionPointer}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.expression.{AnonymousExpressionNode, ExpressionNode, PlainAnonymousExpressionNode, TaskCallInputExpressionNode}
import wom.graph._
import wom.types.{WomOptionalType, WomType}
import wdl.draft3.transforms.wdlom2wdl.WdlWriter.ops._
import wdl.draft3.transforms.wdlom2wdl.WdlWriterImpl.expressionElementWriter

object CallElementToGraphNode {
  def convert(a: CallNodeMakerInputs): ErrorOr[Set[GraphNode]] = {
    val callNodeBuilder = new CallNode.CallNodeBuilder()

    val callName = a.node.alias.getOrElse(a.node.callableReference.split("\\.").last)

    // match the call element to a callable
    def callableValidation: ErrorOr[Callable] =
      a.callables.get(a.node.callableReference) match {
        // pass in specific constructor depending on callable type
        case Some(w: WorkflowDefinition) =>
          val unsuppliedInputs = w.inputs.collect {
            case r: RequiredInputDefinition if r.localName.value.contains(".") => r.localName.value
          }
          val unsuppliedInputsValidation: ErrorOr[Unit] = if (unsuppliedInputs.isEmpty) { ().validNel } else { s"Cannot call '${a.node.callableReference}'. To be called as a sub-workflow it must declare and pass-through the following values via workflow inputs: ${unsuppliedInputs.mkString(", ")}".invalidNel }

          val unspecifiedOutputs = w.graph.outputNodes.map(_.localName).filter(_.contains("."))
          val unspecifiedOutputsValidation: ErrorOr[Unit] = if (unspecifiedOutputs.isEmpty) { ().validNel } else { s"Cannot call '${a.node.callableReference}'. To be called as a sub-workflow it must specify all outputs using an output section. This workflow may wish to declare outputs for: ${unspecifiedOutputs.mkString(", ")}".invalidNel }

          (unsuppliedInputsValidation, unspecifiedOutputsValidation) mapN { (_,_) => w }

        case Some(c: Callable) => c.validNel
        case None => s"Cannot resolve a callable with name ${a.node.callableReference}".invalidNel
      }

    /*
      * Each input definition KV pair becomes an entry in map.
      *
      * i.e.
      * call foo {
      *   input: key = value
      *
      * @return ErrorOr of LocalName(key) mapped to ExpressionNode(value).
      */
    def expressionNodeMappings(callable: Callable): ErrorOr[Map[LocalName, AnonymousExpressionNode]] = {
      def validInput(name: String, definition: Callable.InputDefinition): Boolean = {
        definition.name == name && !definition.isInstanceOf[FixedInputDefinition]
      }

      def hasDeclaration(callable: Callable, name: String): Boolean = callable match {
        case t: TaskDefinition =>
          t.inputs.map(_.name).contains(name)
        case w: WorkflowDefinition =>
          w.graph.nodes.map(_.localName).contains(name)
      }

      a.node.body match {
        case Some(body) =>
          lazy val callNameAlias = a.node.alias match {
            case Some(alias) => s" (as '$alias')"
            case None => ""
          }

          val result = body.inputs.map(input => input.key -> input.value).toMap.traverse { case (name, expression) =>
            callable.inputs.find(i => validInput(name, i)) match {
              case Some(i) =>
                val identifier = WomIdentifier(name)
                val constructor = callable match {
                  case _: CallableTaskDefinition => TaskCallInputExpressionNode.apply _
                  case _ => PlainAnonymousExpressionNode.apply _
                }
                val wdlomWomExpression = WdlomWomExpression(expression, a.linkableValues)
                val requiredInputType = i match {
                  case _: InputDefinitionWithDefault => WomOptionalType(i.womType).flatOptionalType
                  case _ => i.womType
                }

                (WorkflowGraphElementToGraphNode.validateAssignmentType(wdlomWomExpression, requiredInputType) flatMap { _ =>
                  AnonymousExpressionNode.fromInputMapping[AnonymousExpressionNode](identifier, wdlomWomExpression, a.linkablePorts, constructor) map {
                    LocalName(name) -> _
                  }
                }).contextualizeErrors(s"supply input $name = ${expression.toWdlV1}")
              case None =>
                if (hasDeclaration(callable, name)) {
                  s"The call tried to supply a value '$name' that isn't overridable for this task (or sub-workflow). To be able to supply this value, move it into the task (or sub-workflow)'s inputs { } section.".invalidNel
                } else {
                  s"The call supplied a value '$name' that doesn't exist in the task (or sub-workflow)".stripMargin.invalidNel
                }
            }
          }
          result.contextualizeErrors(s"make call to '${callable.name}'$callNameAlias")

        case None => Map.empty[LocalName, AnonymousExpressionNode].valid
      }
    }

    /*
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

        // No input mapping, add an optional input using the default expression
        case withDefault@InputDefinitionWithDefault(n, womType, expression, _, _) =>
          val identifier = WomIdentifier(s"${a.workflowName}.$callName.${n.value}")
          withGraphInputNode(withDefault, OptionalGraphInputNodeWithDefault(identifier, womType, expression, identifier.fullyQualifiedName.value))

        // Not an input, use the default expression:
        case fixedExpression @ FixedInputDefinition(_,_,expression,_, _) => InputDefinitionFold(
          mappings = List(fixedExpression -> Coproduct[InputDefinitionPointer](expression))
        )

        // No input mapping, required and we don't have a default value, create a new RequiredGraphInputNode
        // so that it can be satisfied via workflow inputs
        case required@RequiredInputDefinition(n, womType, _, _) =>
          val identifier = WomIdentifier(s"${a.workflowName}.$callName.${n.value}")
          withGraphInputNode(required, RequiredGraphInputNode(identifier, womType, identifier.fullyQualifiedName.value))

        // No input mapping, no default value but optional, create a OptionalGraphInputNode
        // so that it can be satisfied via workflow inputs
        case optional@OptionalInputDefinition(n, womType, _, _) =>
          val identifier = WomIdentifier(s"${a.workflowName}.$callName.${n.value}")
          withGraphInputNode(optional, OptionalGraphInputNode(identifier, womType, identifier.fullyQualifiedName.value))
      }
    }

    def updateTaskCallNodeInputs(callNodeAndNewNodes: CallNodeAndNewNodes, mappings: Map[LocalName, AnonymousExpressionNode]): Unit = {
      for {
        taskCallNode <- List(callNodeAndNewNodes.node) collect { case c: CommandCallNode => c }
        taskCallInputExpression <- mappings.values.toList collect { case t: TaskCallInputExpressionNode => t }
        _ = taskCallInputExpression.taskCallNodeReceivingInput._graphNode = taskCallNode
      } yield ()
      ()
    }

    for {
      callable <- callableValidation
      mappings <- expressionNodeMappings(callable)
      identifier = WomIdentifier(localName = callName, fullyQualifiedName = a.workflowName + "." + callName)
      result = callNodeBuilder.build(identifier, callable, foldInputDefinitions(mappings, callable))
      _ = updateTaskCallNodeInputs(result, mappings)
    } yield result.nodes
  }
}

case class CallNodeMakerInputs(node: CallElement,
                               linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                               linkablePorts: Map[String, OutputPort],
                               availableTypeAliases: Map[String, WomType],
                               workflowName: String,
                               insideAnotherScatter: Boolean,
                               callables: Map[String, Callable])
