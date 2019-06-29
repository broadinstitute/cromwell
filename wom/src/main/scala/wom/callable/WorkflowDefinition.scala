package wom.callable

import wom.SourceFileLocation
import wom.graph.GraphNode._
import wom.graph.{CommandCallNode, Graph}

final case class WorkflowDefinition(name: String,
                                    innerGraph: Graph,
                                    meta: Map[String, MetaValueElement],
                                    parameterMeta: Map[String, MetaValueElement],
                                    override val sourceLocation : Option[SourceFileLocation]) extends ExecutableCallable {

  override lazy val toString = s"[Workflow $name]"
  override val graph: Graph = innerGraph

  // FIXME: how to get a meaningful order from the node set ?
  override lazy val inputs: List[_ <: Callable.InputDefinition] = innerGraph.nodes.inputDefinitions.toList

  override lazy val outputs: List[_ <: Callable.OutputDefinition] = innerGraph.nodes.outputDefinitions.toList

  override lazy val taskCallNodes: Set[CommandCallNode] = innerGraph.allNodes collect {
    case taskNode: CommandCallNode => taskNode
  }
}
