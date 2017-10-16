package wom.graph

import cats.instances.list._
import cats.kernel.Monoid
import cats.syntax.foldable._
import cats.syntax.traverse._
import lenthall.validation.ErrorOr.ErrorOr
import shapeless.{:+:, CNil, Coproduct}
import wom.callable.Callable._
import wom.callable.{Callable, TaskDefinition, WorkflowDefinition}
import wom.expression.WomExpression
import wom.graph.CallNode._
import wom.graph.GraphNode.GeneratedNodeAndNewNodes
import wom.graph.GraphNodePort.{ConnectedInputPort, GraphNodeOutputPort, InputPort, OutputPort}
import wom.values.WomValue

sealed abstract class CallNode extends GraphNode {
  def callable: Callable
  def callType: String

  def inputDefinitionMappings: InputDefinitionMappings
}

final case class TaskCallNode private(override val identifier: WomIdentifier,
                                      callable: TaskDefinition,
                                      override val inputPorts: Set[GraphNodePort.InputPort],
                                      inputDefinitionMappings: InputDefinitionMappings) extends CallNode {
  val callType: String = "task"
  override val outputPorts: Set[GraphNodePort.OutputPort] = {
    callable.outputs.map(o => GraphNodeOutputPort(o.name, o.womType, this)).toSet
  }
}

final case class WorkflowCallNode private(override val identifier: WomIdentifier,
                                          callable: WorkflowDefinition,
                                          override val inputPorts: Set[GraphNodePort.InputPort],
                                          inputDefinitionMappings: InputDefinitionMappings) extends CallNode {
  val callType: String = "workflow"
  override val outputPorts: Set[GraphNodePort.OutputPort] = {
    callable.innerGraph.nodes.collect { case gon: GraphOutputNode => GraphNodeOutputPort(gon.localName, gon.womType, this) }
  }
}

object TaskCall {
  def graphFromDefinition(taskDefinition: TaskDefinition): ErrorOr[Graph] = {
    val taskDefinitionLocalName = LocalName(taskDefinition.name)
    
    // Creates an identifier for an input or an output
    // The localName is the name of the input or output
    // The FQN combines the name of the task to the name of the input or output
    def identifier(name: LocalName) = WomIdentifier(name, taskDefinitionLocalName.combineToFullyQualifiedName(name))

    def linkOutput(call: GraphNode)(output: OutputDefinition): ErrorOr[GraphNode] = call.outputByName(output.name).map(out => PortBasedGraphOutputNode(
      identifier(output.localName), output.womType, out
    ))
    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap

    val callNodeBuilder = new CallNodeBuilder()
    
    val inputDefinitionFold = taskDefinition.inputs.foldMap({ inputDef =>
    {
      val newNode = inputDef match {
        case RequiredInputDefinition(name, womType) => RequiredGraphInputNode(identifier(name), womType)
        case InputDefinitionWithDefault(name, womType, default) => OptionalGraphInputNodeWithDefault(identifier(name), womType, default)
        case OptionalInputDefinition(name, womType) => OptionalGraphInputNode(identifier(name), womType)
      }

      InputDefinitionFold(
        mappings = Map(inputDef -> Coproduct[InputDefinitionPointer](newNode.singleOutputPort: OutputPort)),
        newGraphInputNodes = Set(newNode),
        callInputPorts = Set(callNodeBuilder.makeInputPort(inputDef, newNode.singleOutputPort))
      )
    }
    })(inputDefinitionFoldMonoid)

    val uniqueIdentifier = WomIdentifier(taskDefinition.name)
    val callWithInputs = callNodeBuilder.build(uniqueIdentifier, taskDefinition, inputDefinitionFold)

    for {
      outputs <- taskDefinition.outputs.traverse(linkOutput(callWithInputs.node) _)
      callSet = Set[GraphNode](callWithInputs.node)
      inputsSet = callWithInputs.newInputs.toSet[GraphNode]
      outputsSet = outputs.toSet[GraphNode]
      graph <- Graph.validateAndConstruct(callSet ++ inputsSet ++ outputsSet)
    } yield graph
  }
}

object CallNode {
  // A monoid can't be derived automatically for this class because it contains a Map[InputDefinition, InputDefinitionPointer],
  // and there's no monoid defined over InputDefinitionPointer
  implicit val inputDefinitionFoldMonoid = new Monoid[InputDefinitionFold] {
    override def empty: InputDefinitionFold = InputDefinitionFold()
    override def combine(x: InputDefinitionFold, y: InputDefinitionFold): InputDefinitionFold = {
      InputDefinitionFold(
        mappings = x.mappings ++ y.mappings,
        callInputPorts = x.callInputPorts ++ y.callInputPorts,
        newGraphInputNodes = x.newGraphInputNodes ++ y.newGraphInputNodes,
        newExpressionNodes = x.newExpressionNodes ++ y.newExpressionNodes
      )
    }
  }

  final case class InputDefinitionFold(mappings: InputDefinitionMappings = Map.empty,
                                                       callInputPorts: Set[InputPort] = Set.empty,
                                                       newGraphInputNodes: Set[GraphInputNode] = Set.empty,
                                                       newExpressionNodes: Set[ExpressionNode] = Set.empty)

  type InputDefinitionPointer = OutputPort :+: WomExpression :+: WomValue :+: CNil
  type InputDefinitionMappings = Map[InputDefinition, InputDefinitionPointer]

  final case class CallNodeAndNewNodes(node: CallNode, newInputs: Set[GraphInputNode], newExpressions: Set[ExpressionNode]) extends GeneratedNodeAndNewNodes {
    def nodes: Set[GraphNode] = Set(node) ++ newInputs ++ newExpressions
  }

  /**
    * Don't use this directly; go via callWithInputs to make sure everything's in order when constructing a CallNode.
    */
  private[graph] def apply(nodeIdentifier: WomIdentifier,
                           callable: Callable,
                           inputPorts: Set[GraphNodePort.InputPort],
                           inputDefinitionMappings: InputDefinitionMappings): CallNode = callable match {
    case t: TaskDefinition => TaskCallNode(nodeIdentifier, t, inputPorts, inputDefinitionMappings)
    case w: WorkflowDefinition => WorkflowCallNode(nodeIdentifier, w, inputPorts, inputDefinitionMappings)
  }

  /**
    * Helper class to build call nodes.
    * Helps making input ports and building the node while making sure node references are set properly.
    */
  class CallNodeBuilder {
    private val graphNodeSetter = new GraphNode.GraphNodeSetter()

    /**
      * Makes an input port for this call.
      * Ensures that the port will contain the reference to the node when it gets created.
      */
    def makeInputPort(inputDefinition: InputDefinition, outputPort: OutputPort) = {
      ConnectedInputPort(inputDefinition.name, inputDefinition.womType, outputPort, graphNodeSetter.get)
    }

    def build(nodeIdentifier: WomIdentifier,
              callable: Callable,
              inputDefinitionFold: InputDefinitionFold): CallNodeAndNewNodes = {
      val callNode = CallNode(nodeIdentifier, callable, inputDefinitionFold.callInputPorts, inputDefinitionFold.mappings)
      graphNodeSetter._graphNode = callNode
      CallNodeAndNewNodes(callNode, inputDefinitionFold.newGraphInputNodes, inputDefinitionFold.newExpressionNodes)
    }
  }
}
