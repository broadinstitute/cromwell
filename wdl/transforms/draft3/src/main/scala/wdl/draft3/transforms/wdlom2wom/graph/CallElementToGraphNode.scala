package wdl.draft3.transforms.wdlom2wom.graph

import cats.syntax.validated._
import cats.syntax.either._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{CallElement, IfElement}
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wom.callable.{CallableTaskDefinition, TaskDefinition}
import wom.{callable, graph}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{CommandCallNode, GraphNode, WomIdentifier}
import wom.types.WomType

object CallElementToGraphNode extends Object {
  def convert(a: CallableNodeMakerInputs): ErrorOr[Set[GraphNode]] = {
    val task = a.tasks.map[CallableTaskDefinition](_.filter(_.name == a.node.callableName).head)
    task map { task =>  Set(CommandCallNode.apply(WomIdentifier(a.node.callableName), task, Set.empty, List.empty)) }
  }
}

final case class CallableNodeMakerInputs(node: CallElement,
                                         tasks: ErrorOr[Set[CallableTaskDefinition]],
                                         linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                         linkablePorts: Map[String, OutputPort],
                                         availableTypeAliases: Map[String, WomType],
                                         workflowName: String,
                                         insideAnotherScatter: Boolean)