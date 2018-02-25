package wdl.draft3.transforms.linking.graph

import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.WorkflowGraphElement
import wdl.model.draft3.graph._
import wdl.model.draft3.graph.UnlinkedValueGenerator.ops._
import wdl.model.draft3.graph.UnlinkedValueConsumer.ops._
import wom.types.WomType

object LinkedGraphMaker {
  def make(nodes: Set[WorkflowGraphElement],
           typeAliases: Map[String, WomType]): ErrorOr[LinkedGraph] = {

    val generatedValuesByGraphNodeValidation = nodes.toList.traverse[ErrorOr, (WorkflowGraphElement, Set[GeneratedValueHandle])] { node =>
      node.generatedValueHandles(typeAliases).map(node -> _)
    } map (_.toMap)

    val consumedValuesByGraphNode = nodes.toList.map(n => n -> n.consumedValueHooks).toMap

    for {
      generatedValuesByGraphNode <- generatedValuesByGraphNodeValidation
      graphNodeByGeneratedValue <- reverseMap(generatedValuesByGraphNode)
      consumedValueLookup <- makeConsumedValueLookup(nodes, graphNodeByGeneratedValue.keySet)
      edges = makeEdges(nodes, consumedValuesByGraphNode, consumedValueLookup, graphNodeByGeneratedValue)
    } yield LinkedGraph(nodes, generatedValuesByGraphNode, edges, consumedValuesByGraphNode, graphNodeByGeneratedValue, consumedValueLookup, typeAliases)
  }

  private def makeEdges(elements: Set[WorkflowGraphElement],
                    consumedValuesByGraphNode: Map[WorkflowGraphElement, Set[UnlinkedConsumedValueHook]],
                    consumedValueLookup: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                    graphElementByGeneratedValueHandle: Map[GeneratedValueHandle, WorkflowGraphElement]): Set[LinkedGraphEdge] = for {
    downstreamElement <- elements
    hook <- consumedValuesByGraphNode(downstreamElement)
    upstreamHandle = consumedValueLookup(hook)
    upstreamElement = graphElementByGeneratedValueHandle(upstreamHandle)
  } yield LinkedGraphEdge(upstreamElement, downstreamElement)

  private def makeConsumedValueLookup(nodes: Set[WorkflowGraphElement],
                                      availableHandles: Set[GeneratedValueHandle]
                                     ): ErrorOr[Map[UnlinkedConsumedValueHook, GeneratedValueHandle]] = {
    val consumed: Set[UnlinkedConsumedValueHook] = nodes.flatMap(_.consumedValueHooks)

    def isMatch(hook: UnlinkedConsumedValueHook, handle: GeneratedValueHandle): Boolean = (hook, handle) match {
      case (UnlinkedIdentifierHook(id1), GeneratedIdentifierValueHandle(id2, _)) => id1 == id2
      case (UnlinkedCallOutputOrIdentifierAndMemberAccessHook(first, _), GeneratedIdentifierValueHandle(id2, _)) if first == id2 => true
      case (UnlinkedCallOutputOrIdentifierAndMemberAccessHook(first1, second1), GeneratedCallOutputValueHandle(first2, second2, _)) if first1 == first2 && second1 == second2 => true
      case _ => false
    }

    def findHandle(consumedValueHook: UnlinkedConsumedValueHook): ErrorOr[(UnlinkedConsumedValueHook, GeneratedValueHandle)] = {
      availableHandles collectFirst {
        case handle if isMatch(consumedValueHook, handle) => handle
      } match {
        case Some(foundHandle) => (consumedValueHook -> foundHandle).validNel
        case None => s"No generated value handle found for consumed value '$consumedValueHook'".invalidNel
      }
    }

    consumed.toList.traverse[ErrorOr, (UnlinkedConsumedValueHook, GeneratedValueHandle)] { findHandle } map {_.toMap}
  }

  private def reverseMap(mapping: Map[WorkflowGraphElement, Set[GeneratedValueHandle]]): ErrorOr[Map[GeneratedValueHandle, WorkflowGraphElement]] = {
    val reversed = for {
      nodeAndHandles <- mapping.toList
      node = nodeAndHandles._1
      handle <- nodeAndHandles._2
    } yield handle -> node

    val reversedMap = reversed.toMap

    if (reversed.lengthCompare(reversedMap.size) == 0) {
      reversedMap.validNel
    } else {
      s"Duplicated generated value handles found in ${reversedMap.keys.map(_.linkableName)}".invalidNel
    }
  }
}
