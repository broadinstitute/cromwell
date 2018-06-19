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
import wom.callable.Callable
import wom.types.WomType

import scalax.collection.Graph
import scalax.collection.GraphEdge.DiEdge

object LinkedGraphMaker {
  def make(nodes: Set[WorkflowGraphElement],
           externalHandles: Set[GeneratedValueHandle],
           typeAliases: Map[String, WomType],
           callables: Map[String, Callable]): ErrorOr[LinkedGraph] = {

    val generatedValuesByGraphNodeValidation = nodes.toList.traverse{ node =>
      node.generatedValueHandles(typeAliases, callables).map(node -> _)
    } map (_.toMap)

    val consumedValuesByGraphNodeValidation: ErrorOr[Map[WorkflowGraphElement, Set[UnlinkedConsumedValueHook]]] = nodes.toList.traverse(n => n.graphElementConsumedValueHooks(typeAliases, callables).map(n -> _)).map(_.toMap)

    for {
      generatedValuesByGraphNode <- generatedValuesByGraphNodeValidation
      consumedValuesByGraphNode <- consumedValuesByGraphNodeValidation
      graphNodeByGeneratedValue <- reverseMap(generatedValuesByGraphNode)
      allHandles = graphNodeByGeneratedValue.keySet ++ externalHandles
      consumedValueLookup <- makeConsumedValueLookup(nodes, typeAliases, allHandles, callables)
      edges = makeEdges(nodes, consumedValuesByGraphNode, consumedValueLookup, graphNodeByGeneratedValue)
    } yield LinkedGraph(nodes, edges, allHandles, consumedValueLookup, typeAliases)
  }

  def getOrdering(linkedGraph: LinkedGraph): ErrorOr[List[WorkflowGraphElement]] = {
    // Find the topological order in which we must create the graph nodes:
    val edges = linkedGraph.edges map { case LinkedGraphEdge(from, to) => DiEdge(from, to) }

    Graph.from[WorkflowGraphElement, DiEdge](linkedGraph.elements, edges).topologicalSort match {
      case Left(cycleNode) => s"This workflow contains a cyclic dependency on ${cycleNode.value}".invalidNel
      // This asInstanceOf is not required, but it suppresses an incorrect intelliJ error highlight:
      case Right(topologicalOrder) => topologicalOrder.toList.map(_.value).asInstanceOf[List[WorkflowGraphElement]].validNel
    }
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
                                      availableHandles: Set[GeneratedValueHandle],
                                      callables: Map[String, Callable]
                                     ): ErrorOr[Map[UnlinkedConsumedValueHook, GeneratedValueHandle]] = {
    val consumedValidation: ErrorOr[Set[UnlinkedConsumedValueHook]] = nodes.toList.traverse(n => n.graphElementConsumedValueHooks(typeAliases, callables)).map(_.toSet.flatten)

    consumedValidation.flatMap { consumed => makeConsumedValueLookup(consumed, availableHandles) }
  }

  def makeConsumedValueLookup(consumedValues: Set[UnlinkedConsumedValueHook], availableHandles: Set[GeneratedValueHandle]): ErrorOr[Map[UnlinkedConsumedValueHook, GeneratedValueHandle]] = {
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
        case None =>
          val didYouMean = availableHandles.map(h => s"'${h.linkableName}'").mkString("[", ", ", "]")
          s"Value '${consumedValueHook.linkString}' is never declared. Available values are: $didYouMean".invalidNel
      }
    }

    consumedValues.toList.traverse { findHandle } map {_.toMap}

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
      s"Duplicate generated values found in: ${reversed.map(_._1.linkableName).sorted.mkString(", ")}".invalidNel
    }
  }
}
