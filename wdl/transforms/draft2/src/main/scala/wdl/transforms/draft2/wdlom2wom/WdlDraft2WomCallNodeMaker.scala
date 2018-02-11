package wdl.transforms.draft2.wdlom2wom

import cats.instances.list._
import cats.syntax.foldable._
import cats.data.NonEmptyList
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import shapeless.Coproduct
import wdl.draft2.model.{WdlCall, WdlTaskCall, WdlWomExpression}
import wdl.draft2.model.{WdlTaskCall, WdlWomExpression}
import wom.callable.Callable
import wom.graph.CallNode._
import wom.callable.Callable.{InputDefinition, InputDefinitionWithDefault, OptionalInputDefinition, RequiredInputDefinition}
import wom.transforms.WomCallableMaker.ops._
import wom.graph.CallNode.{CallNodeAndNewNodes, InputDefinitionFold, InputDefinitionPointer}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.expression.{AnonymousExpressionNode, ExpressionNode, PlainAnonymousExpressionNode, TaskCallInputExpressionNode}
import wom.graph._
import wom.transforms.WomCallNodeMaker

object WdlDraft2WomCallNodeMaker extends WomCallNodeMaker[WdlCall] {
  override def toWomCallNode(taskCall: WdlCall, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, GraphNodePort.OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[CallNode.CallNodeAndNewNodes] = {
    buildWomNodeAndInputs(taskCall, localLookup, outerLookup, preserveIndexForOuterLookups)
  }

  private def buildWomNodeAndInputs(wdlCall: WdlCall, localLookup: Map[String, GraphNodePort.OutputPort], outerLookup: Map[String, GraphNodePort.OutputPort], preserveIndexForOuterLookups: Boolean): ErrorOr[CallNodeAndNewNodes] = {
    import common.validation.ErrorOr._

    val callNodeBuilder = new CallNode.CallNodeBuilder()

    // A validation that all inputs to the call were actually wanted:
    def allInputsWereWantedValidation(callable: Callable): ErrorOr[Unit] = {
      val callableExpectedInputs = callable.inputs.map(_.localName.value)
      val unexpectedInputs: Option[NonEmptyList[String]] = NonEmptyList.fromList(wdlCall.inputMappings.toList collect {
        case (inputName, _) if !callableExpectedInputs.contains(inputName) => inputName
      })

      unexpectedInputs match {
        case None => ().validNel
        case Some(unexpectedInputsNel) => (unexpectedInputsNel map { unexpectedInput =>
          s"Invalid call to '${callable.name}': Didn't expect the input '$unexpectedInput'. Check that this input is declared in the task or workflow. Note that sub-workflow declarations with values that depend on previous values cannot be overridden."
        }).invalid
      }
    }

    /*
      * Each input mapping gets its own ExpressionNode:
      *
      * call my_task { input:
      *   input1 = "hi!"                -> ExpressionNode with no input port
      *   input3 = other_task.out + 2   -> ExpressionNode with an input port pointing to the output port of other_task.out
      * }
     */
    def expressionNodeMappings: ErrorOr[Map[LocalName, AnonymousExpressionNode]] = {
      val precomputedOgins: Map[String, OutputPort] = outerLookup collect {
        case (name, port) if !localLookup.contains(name) => name -> OuterGraphInputNode(WomIdentifier(name), port, preserveIndexForOuterLookups).singleOutputPort
      }
      val newLocalLookup = localLookup ++ precomputedOgins
      wdlCall.inputMappings traverse {
        case (inputName, wdlExpression) =>
          val identifier = wdlCall.womIdentifier.combine(inputName)
          val constructor = wdlCall match {
            case _: WdlTaskCall => TaskCallInputExpressionNode.apply _
            case _ => PlainAnonymousExpressionNode.apply _
          }

          WdlWomExpression.toAnonymousExpressionNode(identifier, WdlWomExpression(wdlExpression, wdlCall), newLocalLookup, Map.empty, preserveIndexForOuterLookups, wdlCall, constructor) map {
            LocalName(inputName) -> _
          }
      }
    }

    /*
      * Fold over the input definitions and
      * 1) assign each input definition its InputDefinitionPointer
      * 2) if necessary, create a graph input node and assign its output port to the input definition
      *
      * The InputDefinitionFold accumulates the input definition mappings, the create graph input nodes, and the expression nodes.
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

      callable.inputs.foldMap {
        // If there is an input mapping for this input definition, use that
        case inputDefinition if expressionNodes.contains(inputDefinition.localName) =>
          val expressionNode = expressionNodes(inputDefinition.localName)
          InputDefinitionFold(
            mappings = List(inputDefinition -> expressionNode.inputDefinitionPointer),
            callInputPorts = Set(callNodeBuilder.makeInputPort(inputDefinition, expressionNode.singleExpressionOutputPort)),
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
          val identifier = wdlCall.womIdentifier.combine(n)
          withGraphInputNode(required, RequiredGraphInputNode(identifier, womType, identifier.fullyQualifiedName.value))

        // No input mapping, no default value but optional, create a OptionalGraphInputNode
        // so that it can be satisfied via workflow inputs
        case optional@OptionalInputDefinition(n, womType, _) =>
          val identifier = wdlCall.womIdentifier.combine(n)
          withGraphInputNode(optional, OptionalGraphInputNode(identifier, womType, identifier.fullyQualifiedName.value))
      }
    }

    (expressionNodeMappings, wdlCall.callable.toWomCallable) flatMapN { case (mappings, callable) =>
      allInputsWereWantedValidation(callable) map { _ =>
        val usedOgins: Set[OuterGraphInputNode] = for {
          expressionNode <- mappings.values.toSet[ExpressionNode]
          ogin <- expressionNode.upstreamOuterGraphInputNodes
        } yield ogin

        val callNodeAndNewNodes = callNodeBuilder.build(wdlCall.womIdentifier, callable, foldInputDefinitions(mappings, callable).copy(usedOuterGraphInputNodes = usedOgins))

        // If the created node is a `TaskCallNode` the created input expressions should be `TaskCallInputExpressionNode`s
        // and should be assigned a reference to the `TaskCallNode`. This is used in the `WorkflowExecutionActor` to
        // find the task and the task's backend so the right `IoFunctionSet` can be used to evaluate task call inputs.
        for {
          taskCallNode <- List(callNodeAndNewNodes.node) collect { case c: CommandCallNode => c }
          taskCallInputExpression <- mappings.values.toList collect { case t: TaskCallInputExpressionNode => t }
          _ = taskCallInputExpression.taskCallNodeReceivingInput._graphNode = taskCallNode
        } yield ()

        callNodeAndNewNodes
      }
    }
  }
}
