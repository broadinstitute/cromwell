package wdl.transforms.base.linking.graph

import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{ExpressionElement, WorkflowGraphElement}
import wdl.model.draft3.graph._
import wdl.model.draft3.graph.UnlinkedValueGenerator.ops._
import wdl.model.draft3.graph.GraphElementValueConsumer.ops._
import wom.callable.Callable
import wom.types.WomType
import scalax.collection.Graph
import scalax.collection.GraphEdge.DiEdge
import wdl.transforms.base.wdlom2wdl.WdlWriter.ops._
import wdl.transforms.base.wdlom2wdl.WdlWriterImpl.graphElementWriter

object LinkedGraphMaker {
  def make(nodes: Set[WorkflowGraphElement],
           externalHandles: Set[GeneratedValueHandle],
           typeAliases: Map[String, WomType],
           callables: Map[String, Callable])
          (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): ErrorOr[LinkedGraph] = {

    val generatedValuesByGraphNodeValidation: ErrorOr[Map[WorkflowGraphElement, Set[GeneratedValueHandle]]] = nodes.toList.traverse{ node =>
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

    def nodeName(workflowGraphElement: WorkflowGraphElement): String = workflowGraphElement.toWdlV1.linesIterator.toList.headOption.getOrElse("Unnamed Element").replace("\"", "")

    // Find the topological order in which we must create the graph nodes:
    val edges = linkedGraph.edges map { case LinkedGraphEdge(from, to) => DiEdge(from, to) }

    val graph = Graph.from[WorkflowGraphElement, DiEdge](linkedGraph.elements, edges)
    graph.topologicalSort match {
      case Left(_) =>
        graph.findCycle match {
          case Some(cycle) if cycle.nonEmpty =>
            // we need to sort the edges lexicographically to make results deterministic. This helps testing.
            // In a cycle like:
            //   b -> d
            //   d -> c
            //   c -> a
            //   a -> b
            // we want to start the cycle with the edge "a -> b"
            val edgeDict : Map[String, String] =
              cycle.value.edges.map{
                case graph.EdgeT(from, to) => nodeName(from) -> nodeName(to)
              }.toMap
            val startPoint = edgeDict.keys.toVector.sorted.head
            var cursor = startPoint

            val cycleReport = scala.collection.mutable.Queue.empty[String]
            // The cycle has the same number of edges as the size of the graph. This gives
            // a clear and simple stopping condition. We do not want to risk
            // a bug that causes an inifinite loop here.
            for (_ <- edgeDict.keys) {
              val next = edgeDict(cursor)
              cycleReport += s""""${cursor}" -> "${next}""""
              cursor = next
            }

            s"""This workflow contains a cyclic dependency:
               |${cycleReport.mkString(System.lineSeparator)}""".stripMargin.invalidNel

          case _ =>
            val edgeStrings = linkedGraph.edges map { case LinkedGraphEdge(from, to) => s""""${nodeName(from)}" -> "${nodeName(to)}"""" }
            // sort the edges for determinism
            val edges = edgeStrings.toVector.sorted
            s"""This workflow contains an elusive cyclic dependency amongst these edges:
               |${edges.mkString(System.lineSeparator)}""".stripMargin.invalidNel
        }

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
                                     )
                                     (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): ErrorOr[Map[UnlinkedConsumedValueHook, GeneratedValueHandle]] = {
    val consumedValidation: ErrorOr[Set[UnlinkedConsumedValueHook]] = nodes.toList.traverse(n => n.graphElementConsumedValueHooks(typeAliases, callables)).map(_.toSet.flatten)

    consumedValidation.flatMap { consumed => makeConsumedValueLookup(consumed, availableHandles) }
  }

  def makeConsumedValueLookup(consumedValues: Set[UnlinkedConsumedValueHook], availableHandles: Set[GeneratedValueHandle]): ErrorOr[Map[UnlinkedConsumedValueHook, GeneratedValueHandle]] = {
    def isMatch(hook: UnlinkedConsumedValueHook, handle: GeneratedValueHandle): Boolean = (hook, handle) match {
      case (UnlinkedIdentifierHook(id1), GeneratedIdentifierValueHandle(id2, _)) => id1 == id2
      case (UnlinkedCallOutputOrIdentifierAndMemberAccessHook(first, _), GeneratedIdentifierValueHandle(id2, _)) if first == id2 => true
      case (UnlinkedCallOutputOrIdentifierAndMemberAccessHook(first1, second1), GeneratedCallOutputValueHandle(first2, second2, _)) if first1 == first2 && second1 == second2 => true
      case (UnlinkedAfterCallHook(upstreamCallName), GeneratedCallFinishedHandle(finishedCallName)) if finishedCallName == upstreamCallName => true
      case _ => false
    }

    def findHandle(consumedValueHook: UnlinkedConsumedValueHook): ErrorOr[(UnlinkedConsumedValueHook, GeneratedValueHandle)] = {
      val maybeFoundHandle = availableHandles collectFirst {
        case handle if isMatch(consumedValueHook, handle) => handle
      }

      (maybeFoundHandle, consumedValueHook) match {
        case (Some(handle), hook) => (hook -> handle).validNel
        case (None, UnlinkedAfterCallHook(upstreamCallName)) =>
            val didYouMean = availableHandles.collect {
              case after: GeneratedCallFinishedHandle => s"'${after.finishedCallName}'"
            }.mkString("[", ", ", "]")
            s"Cannot specify 'after $upstreamCallName': no such call exists. Available calls are: $didYouMean".invalidNel
        case (None, _) =>
            val didYouMean = availableHandles.map(h => s"'${h.linkableName}'").mkString("[", ", ", "]")
            s"Cannot lookup value '${consumedValueHook.linkString}', it is never declared. Available values are: $didYouMean".invalidNel
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
