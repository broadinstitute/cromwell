package wdl.draft3.transforms.wdlom2wom

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.{WorkflowDefinitionElement, WorkflowGraphNodeElement}
import wdl.draft3.transforms.wdlom2wom.ValueConsumer.ops._
import wdl.draft3.transforms.wdlom2wom.ValueGenerator.ops._
import wom.callable.WorkflowDefinition
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{GraphNode, Graph => WomGraph}

import scalax.collection.Graph
import scalax.collection.GraphEdge.DiEdge

object WorkflowDefinitionElementToWomWorkflowDefinition {

  def convert(a: WorkflowDefinitionElement): ErrorOr[WorkflowDefinition] = {

    val graphNodeElements: Set[WorkflowGraphNodeElement] = a.graphElements ++ a.inputsSection.toSeq.flatMap(_.inputDeclarations) ++ a.outputsSection.toSeq.flatMap(_.outputs)

    makeGraph(graphNodeElements, a.name) map { graph =>
      WorkflowDefinition(a.name, graph, Map.empty, Map.empty)
    }

  }

  private def makeGraph(graphNodeElements: Set[WorkflowGraphNodeElement], workflowName: String): ErrorOr[WomGraph] = {
    val unlinkedGraphNodes: Set[UnlinkedGraphNode] = graphNodeElements map UnlinkedGraphNode

    val dependenciesValidation: ErrorOr[List[(UnlinkedGraphNode, UnlinkedGraphNode)]] = generatedByMap(unlinkedGraphNodes) flatMap { gbm =>
      val nodesConsumed = for {
        node <- unlinkedGraphNodes
        consumed <- node.consumedValues
      } yield (node, consumed)

      nodesConsumed.toList.traverse[ErrorOr, (UnlinkedGraphNode, UnlinkedGraphNode)] { case (node, consumed) =>
        findGenerator(consumed, gbm) map { g => (g, node) }
      }
    }

    def graphNodeCreationFold(currentValidation: ErrorOr[List[GraphNode]], next: UnlinkedGraphNode): ErrorOr[List[GraphNode]] = {
      currentValidation flatMap { currentList =>
        val availableValues: Map[String, OutputPort] = (for {
          node <- currentList
          port <- node.outputPorts
        } yield port.name -> port).toMap
        val nextGraphNodeValidation = UnlinkedGraphNodeElementToGraphNode.convert(GraphNodeMakerInputs(next, availableValues, workflowName))
        nextGraphNodeValidation map { nextGraphNode => currentList :+ nextGraphNode }
      }
    }

    val graphNodesValidation = dependenciesValidation flatMap { getOrdering(unlinkedGraphNodes) } flatMap { ordering: List[UnlinkedGraphNode] =>
      ordering.foldLeft[ErrorOr[List[GraphNode]]](List.empty[GraphNode].validNel)(graphNodeCreationFold)
    }

    graphNodesValidation flatMap { graphNodes => WomGraph.validateAndConstruct(graphNodes.toSet) }
  }

  private def getOrdering(unlinkedGraphNodes: Set[UnlinkedGraphNode])
                         (dependencies: List[(UnlinkedGraphNode, UnlinkedGraphNode)]): ErrorOr[List[UnlinkedGraphNode]] = {
    // Find the topological order in which we must create the graph nodes:
    val edges = dependencies map { case (from, to) => DiEdge(from, to) }

    Graph.from[UnlinkedGraphNode, DiEdge](unlinkedGraphNodes, edges).topologicalSort match {
      case Left(cycleNode) => s"This workflow contains a cyclic dependency on ${cycleNode.value}".invalidNel
        // This asInstanceOf is not required, but it suppresses an incorrect intelliJ error highlight:
      case Right(topologicalOrder) => topologicalOrder.toList.map(_.value).asInstanceOf[List[UnlinkedGraphNode]].validNel
    }
  }

  private def findGenerator(consumed: ConsumedValue, generatedByMap: Map[String, UnlinkedGraphNode]): ErrorOr[UnlinkedGraphNode] = consumed match {
    case ConsumedSingleValue(v) if generatedByMap.contains(v) => generatedByMap(v).validNel
    case ConsumedSingleValue(other) => s"Lookup failed for value '$other'".invalidNel
    case ConsumedLookupValue(a, _) if generatedByMap.contains(a) => generatedByMap(a).validNel
    case ConsumedLookupValue(a, b) if generatedByMap.contains(s"$a.$b") => generatedByMap(s"$a.$b").validNel
    case ConsumedLookupValue(a, b) => s"Lookup failed for value '$a.$b'".invalidNel
  }

  private def generatedByMap(unlinkedNodes: Set[UnlinkedGraphNode]): ErrorOr[Map[String, UnlinkedGraphNode]] = {
    def forOneNode(node: UnlinkedGraphNode): ErrorOr[Map[String, UnlinkedGraphNode]] = {
      val generated = node.generatedValues
      val duplicates = generated.groupBy(identity).collect { case (x, l) if l.size > 1 => x }

      val generatedMap = generated.map(x => x -> node).toMap
      if (duplicates.isEmpty) { generatedMap.validNel } else s"Duplicate graph node values created: ${duplicates.mkString(", ")}".invalidNel
    }

    def foldFunction(currentValidation: ErrorOr[Map[String, UnlinkedGraphNode]], node: UnlinkedGraphNode): ErrorOr[Map[String, UnlinkedGraphNode]] = {
      val forThisNodeValidation = forOneNode(node)
      (forThisNodeValidation, currentValidation) flatMapN { (forThisNode, current) =>
        val duplicates = forThisNode.collect { case (k,_) if current.contains(k) => k }
        if (duplicates.isEmpty) {
          (forThisNode ++ current).validNel
        } else {
          s"Duplicate graph node values created: ${duplicates.mkString(", ")}".invalidNel
        }
      }
    }

    val empty: ErrorOr[Map[String, UnlinkedGraphNode]] = Map.empty[String, UnlinkedGraphNode].validNel
    unlinkedNodes.foldLeft[ErrorOr[Map[String, UnlinkedGraphNode]]](empty)(foldFunction)

  }
}
