package wdl.draft3.transforms.wdlom2wom

import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.wdlom2wom.graph.{GraphNodeMakerInputs, WorkflowGraphElementToGraphNode}
import wdl.model.draft3.elements.{WorkflowDefinitionElement, WorkflowGraphElement}
import wdl.draft3.transforms.linking.graph._
import wdl.model.draft3.graph.{GeneratedValueHandle, LinkedGraph, LinkedGraphEdge}
import wom.callable.{Callable, CallableTaskDefinition, TaskDefinition, WorkflowDefinition}
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{CallNode, GraphNode, Graph => WomGraph}
import wom.types.WomType
import scalax.collection.Graph
import scalax.collection.GraphEdge.DiEdge

object WorkflowDefinitionElementToWomWorkflowDefinition {

  final case class WorkflowDefinitionConvertInputs(definitionElement: WorkflowDefinitionElement, typeAliases: Map[String, WomType], callables: Set[Callable])

  def convert(a: WorkflowDefinitionConvertInputs): ErrorOr[WorkflowDefinition] = {

    // Make the set of workflow graph elements, including:
    // - Top-level graph elements
    // - Declarations in the inputs section
    // - Declarations in the outputs section
    val graphNodeElements: Set[WorkflowGraphElement] =
      a.definitionElement.graphElements ++
        a.definitionElement.inputsSection.toSeq.flatMap(_.inputDeclarations) ++
        a.definitionElement.outputsSection.toSeq.flatMap(_.outputs)

    val innerGraph: ErrorOr[WomGraph] = convertGraphElements(GraphLikeConvertInputs(graphNodeElements, Set.empty, a.typeAliases, a.definitionElement.name, insideAScatter = false, a.callables))
    innerGraph map { ig =>  WorkflowDefinition(a.definitionElement.name, ig, Map.empty, Map.empty) }
  }

  final case class GraphLikeConvertInputs(graphElements: Set[WorkflowGraphElement],
                                          seedNodes: Set[GraphNode],
                                          typeAliases: Map[String, WomType],
                                          workflowName: String,
                                          insideAScatter: Boolean,
                                          callables: Set[Callable])

  def convertGraphElements(a: GraphLikeConvertInputs): ErrorOr[WomGraph] = {

    val seedGeneratedValueHandles = for {
      seedNode <- a.seedNodes
      outputPort <- seedNode.outputPorts
    } yield GeneratedValueHandle(outputPort.name, outputPort.womType)

    for {
      linkedGraph <- LinkedGraphMaker.make(nodes = a.graphElements, seedGeneratedValueHandles, typeAliases = a.typeAliases, callables = a.callables)
      womGraph <- makeWomGraph(linkedGraph, a.seedNodes, a.workflowName, a.insideAScatter, a.callables)
    } yield womGraph
  }

  private def makeWomGraph(linkedGraph: LinkedGraph,
                           seedNodes: Set[GraphNode],
                           workflowName: String,
                           insideAScatter: Boolean,
                           callables: Set[Callable]): ErrorOr[WomGraph] = {

    def graphNodeCreationFold(currentValidation: ErrorOr[List[GraphNode]], next: WorkflowGraphElement): ErrorOr[List[GraphNode]] = {
      def outputName(node: GraphNode, port: OutputPort): String = port.identifier.localName.value

      currentValidation flatMap { currentList =>
        val availableValues: Map[String, OutputPort] = (for {
          node <- currentList
          port <- node.outputPorts
        } yield outputName(node, port) -> port).toMap
        val nextGraphNodeValidation = WorkflowGraphElementToGraphNode.convert(GraphNodeMakerInputs(next, linkedGraph.consumedValueLookup, availableValues, linkedGraph.typeAliases, workflowName, insideAScatter, callables))
        nextGraphNodeValidation map { nextGraphNode => currentList ++ nextGraphNode }
      }
    }

    val graphNodesValidation = LinkedGraphMaker.getOrdering(linkedGraph) flatMap { ordering: List[WorkflowGraphElement] =>
      ordering.foldLeft[ErrorOr[List[GraphNode]]](seedNodes.toList.validNel)(graphNodeCreationFold)
    }

    graphNodesValidation flatMap { graphNodes => WomGraph.validateAndConstruct(graphNodes.toSet) }
  }
}
