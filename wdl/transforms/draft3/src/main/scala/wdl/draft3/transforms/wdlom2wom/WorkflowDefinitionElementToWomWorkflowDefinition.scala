package wdl.draft3.transforms.wdlom2wom

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.wdlom2wom.graph.{GraphNodeMakerInputs, WorkflowGraphElementToGraphNode}
import wdl.draft3.transforms.wdlom2wom.linking._
import wdl.model.draft3.elements.{WorkflowDefinitionElement, WorkflowGraphElement}
import wdl.draft3.transforms.wdlom2wom.linking.UnlinkedValueConsumer.ops._
import wdl.draft3.transforms.wdlom2wom.linking.UnlinkedValueGenerator.ops._
import wdl.draft3.transforms.wdlom2wom.graph.linking._
import wom.callable.WorkflowDefinition
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{GraphNode, Graph => WomGraph}

import scalax.collection.Graph
import scalax.collection.GraphEdge.DiEdge

object WorkflowDefinitionElementToWomWorkflowDefinition {

  def convert(a: WorkflowDefinitionElement): ErrorOr[WorkflowDefinition] = {

    val graphNodeElements: Set[WorkflowGraphElement] = a.graphElements ++ a.inputsSection.toSeq.flatMap(_.inputDeclarations) ++ a.outputsSection.toSeq.flatMap(_.outputs)

    makeGraph(graphNodeElements, a.name) map { graph =>
      WorkflowDefinition(a.name, graph, Map.empty, Map.empty)
    }

  }

  private def makeGraph(graphNodeElements: Set[WorkflowGraphElement], workflowName: String): ErrorOr[WomGraph] = {

    val dependenciesValidation: ErrorOr[List[(WorkflowGraphElement, WorkflowGraphElement)]] = generatedByMap(graphNodeElements) flatMap { gbm =>
      val nodesConsumed = for {
        node <- graphNodeElements
        consumed <- node.consumedValueNames
      } yield (node, consumed)

      nodesConsumed.toList.traverse[ErrorOr, (WorkflowGraphElement, WorkflowGraphElement)] { case (node, consumed) =>
        findGenerator(consumed, gbm) map { g => (g, node) }
      }
    }

    def graphNodeCreationFold(currentValidation: ErrorOr[List[GraphNode]], next: WorkflowGraphElement): ErrorOr[List[GraphNode]] = {
      currentValidation flatMap { currentList =>
        val availableValues: Map[String, OutputPort] = (for {
          node <- currentList
          port <- node.outputPorts
        } yield port.name -> port).toMap
        val nextGraphNodeValidation = WorkflowGraphElementToGraphNode.convert(GraphNodeMakerInputs(next, availableValues, workflowName))
        nextGraphNodeValidation map { nextGraphNode => currentList :+ nextGraphNode }
      }
    }

    val graphNodesValidation = dependenciesValidation flatMap { getOrdering(graphNodeElements) } flatMap { ordering: List[WorkflowGraphElement] =>
      ordering.foldLeft[ErrorOr[List[GraphNode]]](List.empty[GraphNode].validNel)(graphNodeCreationFold)
    }

    graphNodesValidation flatMap { graphNodes => WomGraph.validateAndConstruct(graphNodes.toSet) }
  }

  private def getOrdering(unlinkedGraphNodes: Set[WorkflowGraphElement])
                         (dependencies: List[(WorkflowGraphElement, WorkflowGraphElement)]): ErrorOr[List[WorkflowGraphElement]] = {
    // Find the topological order in which we must create the graph nodes:
    val edges = dependencies map { case (from, to) => DiEdge(from, to) }

    Graph.from[WorkflowGraphElement, DiEdge](unlinkedGraphNodes, edges).topologicalSort match {
      case Left(cycleNode) => s"This workflow contains a cyclic dependency on ${cycleNode.value}".invalidNel
        // This asInstanceOf is not required, but it suppresses an incorrect intelliJ error highlight:
      case Right(topologicalOrder) => topologicalOrder.toList.map(_.value).asInstanceOf[List[WorkflowGraphElement]].validNel
    }
  }

  private def findGenerator(consumed: UnlinkedConsumedValueName, generatedByMap: Map[UnlinkedGeneratedValueName, WorkflowGraphElement]): ErrorOr[WorkflowGraphElement] = consumed match {
    case id: UnlinkedIdentifierName if generatedByMap.contains(id) => generatedByMap(id).validNel
    case UnlinkedIdentifierName(other) => s"Lookup failed for value '$other'".invalidNel
    case UnlinkedCallOutputOrIdentifierAndMemberAccess(a, _) if generatedByMap.contains(UnlinkedIdentifierName(a)) =>
      generatedByMap(UnlinkedIdentifierName(a)).validNel
    case UnlinkedCallOutputOrIdentifierAndMemberAccess(a, b) if generatedByMap.contains(UnlinkedCallOutputName(a, b)) =>
      generatedByMap(UnlinkedCallOutputName(a, b)).validNel
    case UnlinkedCallOutputOrIdentifierAndMemberAccess(a, b) => s"Lookup failed for value '$a.$b'".invalidNel
  }

  private def generatedByMap(unlinkedNodes: Set[WorkflowGraphElement]): ErrorOr[Map[UnlinkedGeneratedValueName, WorkflowGraphElement]] = {
    def forOneNode(node: WorkflowGraphElement): ErrorOr[Map[UnlinkedGeneratedValueName, WorkflowGraphElement]] = {
      val generated = node.generatedValueNames
      val duplicates = generated.groupBy(identity).collect { case (x, l) if l.size > 1 => x }

      val generatedMap = generated.map(x => x -> node).toMap
      if (duplicates.isEmpty) { generatedMap.validNel } else s"Duplicate graph node values created: ${duplicates.mkString(", ")}".invalidNel
    }

    def foldFunction(currentValidation: ErrorOr[Map[UnlinkedGeneratedValueName, WorkflowGraphElement]], node: WorkflowGraphElement): ErrorOr[Map[UnlinkedGeneratedValueName, WorkflowGraphElement]] = {
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

    val empty: ErrorOr[Map[UnlinkedGeneratedValueName, WorkflowGraphElement]] = Map.empty[UnlinkedGeneratedValueName, WorkflowGraphElement].validNel
    unlinkedNodes.foldLeft[ErrorOr[Map[UnlinkedGeneratedValueName, WorkflowGraphElement]]](empty)(foldFunction)

  }
}
