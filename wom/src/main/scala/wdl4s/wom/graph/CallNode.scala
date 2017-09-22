package wdl4s.wom.graph

import cats.instances.list._
import cats.kernel.Monoid
import cats.syntax.foldable._
import cats.syntax.traverse._
import lenthall.validation.ErrorOr.ErrorOr
import shapeless.{:+:, CNil, Coproduct}
import wdl4s.wdl.values.WdlValue
import wdl4s.wom.callable.Callable._
import wdl4s.wom.callable.{Callable, TaskDefinition, WorkflowDefinition}
import wdl4s.wom.expression.WomExpression
import wdl4s.wom.graph.CallNode._
import wdl4s.wom.graph.GraphNode.GeneratedNodeAndNewInputs
import wdl4s.wom.graph.GraphNodePort.{ConnectedInputPort, GraphNodeOutputPort, InputPort, OutputPort}

sealed abstract class CallNode extends GraphNode {
  def callable: Callable
  def callType: String

  def inputDefinitionMappings: InputDefinitionMappings
}

final case class TaskCallNode private(override val name: String,
                                      callable: TaskDefinition,
                                      override val inputPorts: Set[GraphNodePort.InputPort],
                                      inputDefinitionMappings: InputDefinitionMappings) extends CallNode {
  val callType: String = "task"
  override val outputPorts: Set[GraphNodePort.OutputPort] = {
    callable.outputs.map(o => GraphNodeOutputPort(o.name, o.womType, this)).toSet
  }
}

final case class WorkflowCallNode private(override val name: String,
                                          callable: WorkflowDefinition,
                                          override val inputPorts: Set[GraphNodePort.InputPort],
                                          inputDefinitionMappings: InputDefinitionMappings) extends CallNode {
  val callType: String = "workflow"
  override val outputPorts: Set[GraphNodePort.OutputPort] = {
    callable.innerGraph.nodes.collect { case gon: GraphOutputNode => GraphNodeOutputPort(gon.name, gon.womType, this) }
  }
}

object TaskCall {
  def graphFromDefinition(taskDefinition: TaskDefinition): ErrorOr[Graph] = {
    def linkOutput(call: GraphNode)(output: OutputDefinition): ErrorOr[GraphNode] = call.outputByName(output.name).map(out => PortBasedGraphOutputNode(output.name, output.womType, out))
    import lenthall.validation.ErrorOr.ShortCircuitingFlatMap

    val callNodeBuilder = new CallNodeBuilder()

    val inputDefinitionFold = taskDefinition.inputs.foldMap({ inputDef =>
    {
      val newNode = inputDef match {
        case RequiredInputDefinition(name, womType) => RequiredGraphInputNode(name, womType)
        case InputDefinitionWithDefault(name, womType, default) => OptionalGraphInputNodeWithDefault(name, womType, default)
        case OptionalInputDefinition(name, womType) => OptionalGraphInputNode(name, womType)
      }

      InputDefinitionFold(
        mappings = Map(inputDef -> Coproduct[InputDefinitionPointer](newNode.singleOutputPort: OutputPort)),
        newGraphInputNodes = Set(newNode),
        callInputPorts = Set(callNodeBuilder.makeInputPort(inputDef, newNode.singleOutputPort))
      )
    }
    })(inputDefinitionFoldMonoid)

    val callWithInputs = callNodeBuilder.build(taskDefinition.name, taskDefinition, inputDefinitionFold)

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

  private [wdl4s] final case class InputDefinitionFold(mappings: InputDefinitionMappings = Map.empty,
                                                       callInputPorts: Set[InputPort] = Set.empty,
                                                       newGraphInputNodes: Set[GraphInputNode] = Set.empty,
                                                       newExpressionNodes: Set[ExpressionNode] = Set.empty)

  type InputDefinitionPointer = OutputPort :+: WomExpression :+: WdlValue :+: CNil
  type InputDefinitionMappings = Map[InputDefinition, InputDefinitionPointer]

  final case class CallNodeAndNewNodes(node: CallNode, newInputs: Set[GraphInputNode], newExpressions: Set[ExpressionNode]) extends GeneratedNodeAndNewInputs {
    def nodes: Set[GraphNode] = Set(node) ++ newInputs ++ newExpressions
  }

  /**
    * Don't use this directly; go via callWithInputs to make sure everything's in order when constructing a CallNode.
    */
  private[graph] def apply(name: String,
                           callable: Callable,
                           inputPorts: Set[GraphNodePort.InputPort],
                           inputDefinitionMappings: InputDefinitionMappings): CallNode = callable match {
    case t: TaskDefinition => TaskCallNode(name, t, inputPorts, inputDefinitionMappings)
    case w: WorkflowDefinition => WorkflowCallNode(name, w, inputPorts, inputDefinitionMappings)
  }

  /**
    * Helper class to build call nodes.
    * Helps making input ports and building the node while making sure node references are set properly.
    */
  private [wdl4s] class CallNodeBuilder {
    private val graphNodeSetter = new GraphNode.GraphNodeSetter()

    /**
      * Makes an input port for this call.
      * Ensures that the port will contain the reference to the node when it gets created.
      */
    def makeInputPort(inputDefinition: InputDefinition, outputPort: OutputPort) = {
      ConnectedInputPort(inputDefinition.name, inputDefinition.womType, outputPort, graphNodeSetter.get)
    }

    def build(name: String,
              callable: Callable,
              inputDefinitionFold: InputDefinitionFold): CallNodeAndNewNodes = {
      val callNode = CallNode(name, callable, inputDefinitionFold.callInputPorts, inputDefinitionFold.mappings)
      graphNodeSetter._graphNode = callNode
      CallNodeAndNewNodes(callNode, inputDefinitionFold.newGraphInputNodes, inputDefinitionFold.newExpressionNodes)
    }
  }
}
