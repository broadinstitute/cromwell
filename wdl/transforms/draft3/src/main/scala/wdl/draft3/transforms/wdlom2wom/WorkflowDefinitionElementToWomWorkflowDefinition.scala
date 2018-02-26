package wdl.draft3.transforms.wdlom2wom

import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.wdlom2wom.graph.{GraphNodeMakerInputs, WorkflowGraphElementToGraphNode}
import wdl.model.draft3.elements.{WorkflowDefinitionElement, WorkflowGraphElement}
import wdl.draft3.transforms.linking.graph._
import wdl.model.draft3.graph.{LinkedGraph, LinkedGraphEdge}
import wom.callable.WorkflowDefinition
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{GraphNode, Graph => WomGraph}

import scalax.collection.Graph
import scalax.collection.GraphEdge.DiEdge

object WorkflowDefinitionElementToWomWorkflowDefinition {

  def convert(a: WorkflowDefinitionElement): ErrorOr[WorkflowDefinition] = {

    val graphNodeElements: Set[WorkflowGraphElement] = a.graphElements ++ a.inputsSection.toSeq.flatMap(_.inputDeclarations) ++ a.outputsSection.toSeq.flatMap(_.outputs)

    for {
      linkedGraph <- LinkedGraphMaker.make(nodes = graphNodeElements, typeAliases = Map.empty)
      womGraph <- makeWomGraph(linkedGraph, a.name)
    } yield WorkflowDefinition(a.name, womGraph, Map.empty, Map.empty)

  }

  private def makeWomGraph(linkedGraph: LinkedGraph, workflowName: String): ErrorOr[WomGraph] = {

    def graphNodeCreationFold(currentValidation: ErrorOr[List[GraphNode]], next: WorkflowGraphElement): ErrorOr[List[GraphNode]] = {
      currentValidation flatMap { currentList =>
        val availableValues: Map[String, OutputPort] = (for {
          node <- currentList
          port <- node.outputPorts
        } yield port.name -> port).toMap
        val nextGraphNodeValidation = WorkflowGraphElementToGraphNode.convert(GraphNodeMakerInputs(next, linkedGraph.consumedValueLookup, availableValues, workflowName))
        nextGraphNodeValidation map { nextGraphNode => currentList :+ nextGraphNode }
      }
    }

    val graphNodesValidation = getOrdering(linkedGraph) flatMap { ordering: List[WorkflowGraphElement] =>
      ordering.foldLeft[ErrorOr[List[GraphNode]]](List.empty[GraphNode].validNel)(graphNodeCreationFold)
    }

    graphNodesValidation flatMap { graphNodes => WomGraph.validateAndConstruct(graphNodes.toSet) }
  }

  private def getOrdering(linkedGraph: LinkedGraph): ErrorOr[List[WorkflowGraphElement]] = {
    // Find the topological order in which we must create the graph nodes:
    val edges = linkedGraph.edges map { case LinkedGraphEdge(from, to) => DiEdge(from, to) }

    Graph.from[WorkflowGraphElement, DiEdge](linkedGraph.elements, edges).topologicalSort match {
      case Left(cycleNode) => s"This workflow contains a cyclic dependency on ${cycleNode.value}".invalidNel
        // This asInstanceOf is not required, but it suppresses an incorrect intelliJ error highlight:
      case Right(topologicalOrder) => topologicalOrder.toList.map(_.value).asInstanceOf[List[WorkflowGraphElement]].validNel
    }
  }
}
