package wom.graph

import cats.data.Validated.Valid
import cats.instances.list._
import cats.kernel.Monoid
import cats.syntax.validated._
import cats.syntax.foldable._
import common.validation.ErrorOr.ErrorOr
import shapeless.{:+:, CNil, Coproduct}
import wom.callable.Callable._
import wom.callable.CommandTaskDefinition.OutputFunctionResponse
import wom.callable._
import wom.expression.{InputPointerToWomValue, IoFunctionSet, WomExpression}
import wom.graph.CallNode._
import wom.graph.GraphNode.GeneratedNodeAndNewNodes
import wom.graph.GraphNodePort._
import wom.graph.expression.{ExpressionNode, ExpressionNodeLike}
import wom.values.{WomEvaluatedCallInputs, WomValue}

import scala.concurrent.ExecutionContext

sealed abstract class CallNode extends GraphNode {
  def callable: Callable
  def callType: String

  def inputDefinitionMappings: InputDefinitionMappings
}

final case class ExpressionCallNode private(override val identifier: WomIdentifier,
                                      callable: ExpressionTaskDefinition,
                                      override val inputPorts: Set[GraphNodePort.InputPort],
                                      inputDefinitionMappings: InputDefinitionMappings,
                                      outputIdentifierCompoundingFunction: (WomIdentifier, String) => WomIdentifier) extends CallNode with ExpressionNodeLike {
  val callType: String = "expression task"
  lazy val expressionBasedOutputPorts: List[ExpressionBasedOutputPort] = {
    callable.outputs.map(o => ExpressionBasedOutputPort(outputIdentifierCompoundingFunction(identifier, o.localName.value), o.womType, this, o.expression))
  }

  override lazy val outputPorts: Set[OutputPort] = expressionBasedOutputPorts.toSet[OutputPort]

  override def evaluate(outputPortLookup: OutputPort => ErrorOr[WomValue], ioFunctionSet: IoFunctionSet) = {
    for {
      // Evaluate the inputs to get a lookup to evaluate the actual expression
      womEvaluatedInputs <- CallNode.resolveAndEvaluateInputs(this, ioFunctionSet, outputPortLookup).toEither
      // Make a usable lookup
      lookup = womEvaluatedInputs.map({ case (inputDefinition, value) => inputDefinition.name -> value })
      // Evaluate the expression
      evaluated <- callable.evaluate(lookup, ioFunctionSet, expressionBasedOutputPorts)
    } yield evaluated
  }
}

final case class CommandCallNode private(override val identifier: WomIdentifier,
                                         callable: CommandTaskDefinition,
                                         override val inputPorts: Set[GraphNodePort.InputPort],
                                         inputDefinitionMappings: InputDefinitionMappings,
                                         outputIdentifierCompoundingFunction: (WomIdentifier, String) => WomIdentifier) extends CallNode {
  val callType: String = "task"
  lazy val expressionBasedOutputPorts: List[ExpressionBasedOutputPort] = {
    callable.outputs.map(o => ExpressionBasedOutputPort(outputIdentifierCompoundingFunction(identifier, o.localName.value), o.womType, this, o.expression))
  }

  override lazy val outputPorts: Set[OutputPort] = expressionBasedOutputPorts.toSet[OutputPort]

  /**
    * Evaluate outputs using the custom evaluation function of the task definition.
    * An empty return value means the engine should fall back to its default evaluation method.
    */
  def customOutputEvaluation(inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, executionContext: ExecutionContext): OutputFunctionResponse = {
    callable.customizedOutputEvaluation(outputPorts, inputs, ioFunctionSet, executionContext)
  }
}

final case class WorkflowCallNode private (override val identifier: WomIdentifier,
                                           callable: WorkflowDefinition,
                                           override val inputPorts: Set[GraphNodePort.InputPort],
                                           inputDefinitionMappings: InputDefinitionMappings,
                                           outputIdentifierCompoundingFunction: (WomIdentifier, String) => WomIdentifier) extends CallNode {
  val callType: String = "workflow"
  val subworkflowCallOutputPorts: Set[SubworkflowCallOutputPort] = {
    callable.innerGraph.nodes.collect { case gon: GraphOutputNode => SubworkflowCallOutputPort(outputIdentifierCompoundingFunction(identifier, gon.localName), gon, this) }
  }
  override val outputPorts: Set[OutputPort] = subworkflowCallOutputPorts.toSet[OutputPort]
}



object TaskCall {
  def graphFromDefinition(taskDefinition: TaskDefinition): ErrorOr[Graph] = {
    val taskDefinitionLocalName = LocalName(taskDefinition.name)
    
    /* Creates an identifier for an input or an output
     * The localName is the name of the input or output
     */
    
    // The FQN combines the name of the task to the name of the input or output
    def identifier(name: LocalName) = WomIdentifier(name, taskDefinitionLocalName.combineToFullyQualifiedName(name))

    import common.validation.ErrorOr.ShortCircuitingFlatMap

    val callNodeBuilder = new CallNodeBuilder()
    
    val inputDefinitionFold = taskDefinition.inputs.foldMap({ inputDef =>
    {
      val newNode: Option[ExternalGraphInputNode] = inputDef match {
        case RequiredInputDefinition(name, womType, _, _) => Some(RequiredGraphInputNode(identifier(name), womType, name.value))
        case InputDefinitionWithDefault(name, womType, default, _, _) => Some(OptionalGraphInputNodeWithDefault(identifier(name), womType, default, name.value))
        case OptionalInputDefinition(name, womType, _, _) => Some(OptionalGraphInputNode(identifier(name), womType, name.value))
        case _: FixedInputDefinition => None
      }

      newNode match {
        case Some(inputNode) =>
          InputDefinitionFold(
            mappings = List(inputDef -> Coproduct[InputDefinitionPointer](inputNode.singleOutputPort: OutputPort)),
            newGraphInputNodes = Set(inputNode),
            callInputPorts = Set(callNodeBuilder.makeInputPort(inputDef, inputNode.singleOutputPort))
          )
        case None => InputDefinitionFold() // No-op
      }

    }
    })(inputDefinitionFoldMonoid)

    val uniqueIdentifier = WomIdentifier(taskDefinition.name)
    val callWithInputs = callNodeBuilder.build(uniqueIdentifier, taskDefinition, inputDefinitionFold)

    val outputNodes: ErrorOr[Seq[GraphOutputNode]] = callWithInputs.node.outputPorts.map(output => PortBasedGraphOutputNode(
      identifier(LocalName(output.internalName)), output.womType, output
    )).toList.validNel
    val result = for {
      outputs <- outputNodes
      callSet = Set[GraphNode](callWithInputs.node)
      inputsSet = callWithInputs.newInputs.toSet[GraphNode]
      outputsSet = outputs.toSet[GraphNode]
      graph <- Graph.validateAndConstruct(callSet ++ inputsSet ++ outputsSet)
    } yield graph

    result
  }
}

object CallNode {
  def resolveAndEvaluateInputs(callNode: CallNode,
                               expressionLanguageFunctions: IoFunctionSet,
                               outputStoreLookup: OutputPort => ErrorOr[WomValue]): ErrorOr[WomEvaluatedCallInputs] = {
    import common.validation.Validation._
    import common.validation.ErrorOr._
    
    callNode.inputDefinitionMappings.foldLeft(Map.empty[InputDefinition, ErrorOr[WomValue]]) {
      case (accumulatedInputsSoFar, (inputDefinition, pointer)) =>
        // We could have a commons method for this kind of "filtering valid values"
        val validInputsAccumulated: Map[String, WomValue] = accumulatedInputsSoFar.collect({
          case (input, Valid(errorOrWdlValue)) => input.name -> errorOrWdlValue
        })

        val coercedValue = pointer.fold(InputPointerToWomValue).apply(
          validInputsAccumulated, expressionLanguageFunctions, outputStoreLookup, inputDefinition
        ) flatMap(inputDefinition.womType.coerceRawValue(_).toErrorOr)

        accumulatedInputsSoFar + (inputDefinition -> coercedValue)
    }.sequence
  }
  
  /* A monoid can't be derived automatically for this class because it contains a Map[InputDefinition, InputDefinitionPointer],
   * and there's no monoid defined over InputDefinitionPointer
   */
  implicit val inputDefinitionFoldMonoid = new Monoid[InputDefinitionFold] {
    override def empty: InputDefinitionFold = InputDefinitionFold()
    override def combine(x: InputDefinitionFold, y: InputDefinitionFold): InputDefinitionFold = {
      InputDefinitionFold(
        mappings = x.mappings ++ y.mappings,
        callInputPorts = x.callInputPorts ++ y.callInputPorts,
        newGraphInputNodes = x.newGraphInputNodes ++ y.newGraphInputNodes,
        newExpressionNodes = x.newExpressionNodes ++ y.newExpressionNodes,
        usedOuterGraphInputNodes = x.usedOuterGraphInputNodes ++ y.usedOuterGraphInputNodes
      )
    }
  }

  final case class InputDefinitionFold(mappings: InputDefinitionMappings = List.empty,
                                       callInputPorts: Set[InputPort] = Set.empty,
                                       newGraphInputNodes: Set[ExternalGraphInputNode] = Set.empty,
                                       newExpressionNodes: Set[ExpressionNode] = Set.empty,
                                       usedOuterGraphInputNodes: Set[OuterGraphInputNode] = Set.empty)

  type InputDefinitionPointer = OutputPort :+: WomExpression :+: WomValue :+: CNil
  // This is a List rather than Map because the order of 'InputDefinition's is important:
  type InputDefinitionMappings = List[(InputDefinition, InputDefinitionPointer)]

  final case class CallNodeAndNewNodes(node: CallNode, newInputs: Set[ExternalGraphInputNode], newExpressions: Set[ExpressionNode], override val usedOuterGraphInputNodes: Set[OuterGraphInputNode]) extends GeneratedNodeAndNewNodes

  /**
    * Don't use this directly; go via callWithInputs to make sure everything's in order when constructing a CallNode.
    */
  private[graph] def apply(nodeIdentifier: WomIdentifier,
                           callable: Callable,
                           inputPorts: Set[GraphNodePort.InputPort],
                           inputDefinitionMappings: InputDefinitionMappings,
                           outputIdentifierCompoundingFunction: (WomIdentifier, String) => WomIdentifier): CallNode = callable match {
    case t: CommandTaskDefinition => CommandCallNode(nodeIdentifier, t, inputPorts, inputDefinitionMappings, outputIdentifierCompoundingFunction)
    case w: WorkflowDefinition => WorkflowCallNode(nodeIdentifier, w, inputPorts, inputDefinitionMappings, outputIdentifierCompoundingFunction)
    case w: ExpressionTaskDefinition => ExpressionCallNode(nodeIdentifier, w, inputPorts, inputDefinitionMappings, outputIdentifierCompoundingFunction)
  }

  /**
    * Helper class to build call nodes.
    * Helps making input ports and building the node while making sure node references are set properly.
    */
  class CallNodeBuilder {
    private val graphNodeSetter = new GraphNode.GraphNodeSetter[CallNode]()

    /**
      * Makes an input port for this call.
      * Ensures that the port will contain the reference to the node when it gets created.
      */
    def makeInputPort(inputDefinition: InputDefinition, outputPort: OutputPort) = {
      ConnectedInputPort(inputDefinition.name, inputDefinition.womType, outputPort, graphNodeSetter.get)
    }

    def build(nodeIdentifier: WomIdentifier,
              callable: Callable,
              inputDefinitionFold: InputDefinitionFold,
              outputIdentifierCompoundingFunction: (WomIdentifier, String) => WomIdentifier = defaultOutputIdentifierCompounder
             ): CallNodeAndNewNodes = {
      val callNode = CallNode(nodeIdentifier, callable, inputDefinitionFold.callInputPorts, inputDefinitionFold.mappings, outputIdentifierCompoundingFunction)
      graphNodeSetter._graphNode = callNode
      CallNodeAndNewNodes(callNode, inputDefinitionFold.newGraphInputNodes, inputDefinitionFold.newExpressionNodes, inputDefinitionFold.usedOuterGraphInputNodes)
    }

    def defaultOutputIdentifierCompounder(callIdentifier: WomIdentifier, outputName: String): WomIdentifier = {
      callIdentifier.combine(outputName)
    }
  }
}
