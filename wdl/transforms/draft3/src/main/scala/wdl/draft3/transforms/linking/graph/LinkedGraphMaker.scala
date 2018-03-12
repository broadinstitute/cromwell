package wdl.draft3.transforms.linking.graph

import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.WorkflowGraphElement
import wdl.model.draft3.graph._
import wdl.model.draft3.graph.UnlinkedValueGenerator.ops._
import wdl.model.draft3.graph.GraphElementValueConsumer.ops._
import wom.types.WomType

object LinkedGraphMaker {
  def make(nodes: Set[WorkflowGraphElement],
           externalHandles: Set[GeneratedValueHandle],
           typeAliases: Map[String, WomType]): ErrorOr[LinkedGraph] = {

    val generatedValuesByGraphNodeValidation = nodes.toList.traverse[ErrorOr, (WorkflowGraphElement, Set[GeneratedValueHandle])] { node =>
      node.generatedValueHandles(typeAliases).map(node -> _)
    } map (_.toMap)

    val consumedValuesByGraphNodeValidation: ErrorOr[Map[WorkflowGraphElement, Set[UnlinkedConsumedValueHook]]] = nodes.toList.traverse[ErrorOr, (WorkflowGraphElement, Set[UnlinkedConsumedValueHook])](n => n.graphElementConsumedValueHooks(typeAliases).map(n -> _)).map(_.toMap)

    for {
      generatedValuesByGraphNode <- generatedValuesByGraphNodeValidation
      consumedValuesByGraphNode <- consumedValuesByGraphNodeValidation
      graphNodeByGeneratedValue <- reverseMap(generatedValuesByGraphNode)
      consumedValueLookup <- makeConsumedValueLookup(nodes, typeAliases, graphNodeByGeneratedValue.keySet ++ externalHandles)
      edges = makeEdges(nodes, consumedValuesByGraphNode, consumedValueLookup, graphNodeByGeneratedValue)
    } yield LinkedGraph(nodes, edges, consumedValueLookup, typeAliases)
  }

  private def makeEdges(elements: Set[WorkflowGraphElement],
                    consumedValuesByGraphNode: Map[WorkflowGraphElement, Set[UnlinkedConsumedValueHook]],
                    consumedValueLookup: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                    graphElementByGeneratedValueHandle: Map[GeneratedValueHandle, WorkflowGraphElement]): Set[LinkedGraphEdge] = for {
    downstreamElement <- elements
    hook <- consumedValuesByGraphNode(downstreamElement)
    upstreamHandle = consumedValueLookup(hook)
    upstreamElement <- graphElementByGeneratedValueHandle.get(upstreamHandle).toList
  } yield LinkedGraphEdge(upstreamElement, downstreamElement)

  private def makeConsumedValueLookup(nodes: Set[WorkflowGraphElement],
                                      typeAliases: Map[String, WomType],
                                      availableHandles: Set[GeneratedValueHandle]
                                     ): ErrorOr[Map[UnlinkedConsumedValueHook, GeneratedValueHandle]] = {
    val consumedValidation: ErrorOr[Set[UnlinkedConsumedValueHook]] = nodes.toList.traverse(n => n.graphElementConsumedValueHooks(typeAliases)).map(_.toSet.flatten)

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

    consumedValidation.flatMap { consumed =>
      consumed.toList.traverse[ErrorOr, (UnlinkedConsumedValueHook, GeneratedValueHandle)] { findHandle } map {_.toMap}
    }
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
