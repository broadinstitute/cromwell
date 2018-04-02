package wdl.draft3.transforms.wdlom2wom.graph

import cats.syntax.validated._
import cats.syntax.either._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{CallElement, IfElement}
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.callable.{Callable, CallableTaskDefinition, CommandTaskDefinition, TaskDefinition}
import wom.graph.CallNode.{CallNodeBuilder, InputDefinitionFold}
import wom.{callable, graph}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{CommandCallNode, GraphNode, WomIdentifier}
import wom.types.WomType

object CallElementToGraphNode extends Object {
  def convert(a: CallableNodeMakerInputs): ErrorOr[Set[GraphNode]] = {
    val callable: ErrorOr[Callable] = a.callables.find(_.name == a.node.callableName) match {
      case Some(task: CommandTaskDefinition) => task.valid
      case Some(c: Callable) => c.valid
      case None => s"Cannot resolve a callable with name ${a.node.callableName}".invalidNel
    }

    callable map { callable =>
      new CallNodeBuilder().build(WomIdentifier(a.node.callableName), callable, InputDefinitionFold(), true).nodes
    }
  }
}

final case class CallableNodeMakerInputs(node: CallElement,
                                         linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                         linkablePorts: Map[String, Any],
                                         availableTypeAliases: Map[String, WomType],
                                         workflowName: String,
                                         insideAnotherScatter: Boolean,
                                         callables: Set[Callable])