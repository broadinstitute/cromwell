package wdl.transforms.base.wdlom2wom

import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, _}
import wdl.model.draft3.elements.ExpressionElement.{ArrayLiteral, IdentifierLookup, SelectFirst}
import wdl.model.draft3.elements._
import wdl.model.draft3.graph.{ExpressionValueConsumer, GeneratedValueHandle, LinkedGraph}
import wdl.shared.transforms.wdlom2wom.WomGraphMakerTools
import wdl.transforms.base.linking.graph.LinkedGraphMaker
import wdl.transforms.base.wdlom2wom.graph.{GraphNodeMakerInputs, WorkflowGraphElementToGraphNode}
import wom.callable.{Callable, WorkflowDefinition}
import wom.graph.expression.AnonymousExpressionNode
import wom.graph.GraphNodePort.OutputPort
import wom.graph.{CallNode, GraphNode, WomIdentifier, Graph => WomGraph}
import wom.types.WomType
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.model.draft3.graph.expression.{FileEvaluator, TypeEvaluator, ValueEvaluator}
import wdl.transforms.base.wdlom2wom.graph.renaming.GraphIdentifierLookupRenamer.ops._
import wdl.transforms.base.wdlom2wom.graph.renaming._

object WorkflowDefinitionElementToWomWorkflowDefinition {

  final case class WorkflowDefinitionConvertInputs(definitionElement: WorkflowDefinitionElement, typeAliases: Map[String, WomType], callables: Map[String, Callable])

  def convert(b: WorkflowDefinitionConvertInputs)
             (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
              fileEvaluator: FileEvaluator[ExpressionElement],
              typeEvaluator: TypeEvaluator[ExpressionElement],
              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[WorkflowDefinition] = {

    val a: WorkflowDefinitionConvertInputs = eliminateInputDependencies(b)

    // Make the set of workflow graph elements, including:
    // - Top-level graph elements
    // - Declarations in the inputs section
    // - Declarations in the outputs section
    val graphNodeElements: Set[WorkflowGraphElement] =
      a.definitionElement.graphElements ++
        a.definitionElement.inputsSection.toSeq.flatMap(_.inputDeclarations) ++
        a.definitionElement.outputsSection.toSeq.flatMap(_.outputs)

    val innerGraph: ErrorOr[WomGraph] = convertGraphElements(GraphLikeConvertInputs(graphNodeElements, Set.empty, a.typeAliases, a.definitionElement.name, insideAScatter = false, a.callables))
    // NB: isEmpty means "not isDefined". We specifically do NOT add defaults if the output section is defined but empty.
    val withDefaultOutputs: ErrorOr[WomGraph] = if (a.definitionElement.outputsSection.isEmpty) {
      innerGraph map { WomGraphMakerTools.addDefaultOutputs(_, Some(WomIdentifier(a.definitionElement.name))) }
    } else {
      innerGraph
    }
    (withDefaultOutputs map { ig =>  WorkflowDefinition(a.definitionElement.name, ig, Map.empty, Map.empty) }).contextualizeErrors(s"process workflow definition '${a.definitionElement.name}'")
  }

  final case class GraphLikeConvertInputs(graphElements: Set[WorkflowGraphElement],
                                          seedNodes: Set[GraphNode],
                                          typeAliases: Map[String, WomType],
                                          workflowName: String,
                                          insideAScatter: Boolean,
                                          callables: Map[String, Callable])

  def convertGraphElements(a: GraphLikeConvertInputs)
                          (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
                           fileEvaluator: FileEvaluator[ExpressionElement],
                           typeEvaluator: TypeEvaluator[ExpressionElement],
                           valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[WomGraph] = {

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
                           callables: Map[String, Callable])
                          (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
                           fileEvaluator: FileEvaluator[ExpressionElement],
                           typeEvaluator: TypeEvaluator[ExpressionElement],
                           valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[WomGraph] = {

    def graphNodeCreationFold(currentValidation: ErrorOr[List[GraphNode]], next: WorkflowGraphElement): ErrorOr[List[GraphNode]] = {
      def outputName(node: GraphNode, port: OutputPort): String = port.identifier.localName.value

      currentValidation flatMap { currentList =>
        val availableValues: Map[String, OutputPort] = (for {
          // Anonymous expression nodes are one-time-use and do not exist in the universe of available linking candidates (#3999)
          node <- currentList.filterNot(_.isInstanceOf[AnonymousExpressionNode])
          port <- node.outputPorts
        } yield outputName(node, port) -> port).toMap

        val upstreamGraphNodes = (currentList collect {
          case node: CallNode => node.localName -> node
        }).toMap

        val generatedGraphNodesValidation: ErrorOr[Set[GraphNode]] =
          WorkflowGraphElementToGraphNode.convert(
            GraphNodeMakerInputs(next, upstreamGraphNodes, linkedGraph.consumedValueLookup, availableValues, linkedGraph.typeAliases, workflowName, insideAScatter, callables))
        generatedGraphNodesValidation map { nextGraphNodes: Set[GraphNode] => currentList ++ nextGraphNodes }
      }
    }

    val graphNodesValidation = LinkedGraphMaker.getOrdering(linkedGraph) flatMap { ordering: List[WorkflowGraphElement] =>
      ordering.foldLeft[ErrorOr[List[GraphNode]]](seedNodes.toList.validNel)(graphNodeCreationFold)
    }

    graphNodesValidation flatMap { graphNodes => WomGraph.validateAndConstruct(graphNodes.toSet) }
  }

  private def eliminateInputDependencies(a: WorkflowDefinitionConvertInputs)
                                        (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement]): WorkflowDefinitionConvertInputs = {
    case class NewInputElementsSet(original: InputDeclarationElement, newInput: InputDeclarationElement, newDeclaration: IntermediateValueDeclarationElement)

    val inputElementsWithUpstreams: Seq[NewInputElementsSet] = a.definitionElement.inputsSection.map(_.inputDeclarations).getOrElse(Seq.empty) collect {
      case ide @ InputDeclarationElement(typeElement,name, Some(expression)) if expression.expressionConsumedValueHooks.nonEmpty =>
        val input = InputDeclarationElement(OptionalTypeElement(typeElement), name, None)

        val selecterExpression = SelectFirst(ArrayLiteral(Seq(IdentifierLookup(name), expression)))
        val intermediate = IntermediateValueDeclarationElement(typeElement, s"__$name", selecterExpression)

        NewInputElementsSet(ide, input, intermediate)
    }

    if (inputElementsWithUpstreams.nonEmpty) {
      val newInputsSection: Option[InputsSectionElement] = a.definitionElement.inputsSection.map { inputsSection =>
        InputsSectionElement(inputsSection.inputDeclarations.filterNot(inputElementsWithUpstreams.map(_.original).contains) ++ inputElementsWithUpstreams.map(_.newInput))
      }
      val newGraphNodeSet: Set[WorkflowGraphElement] = a.definitionElement.graphElements ++ inputElementsWithUpstreams.map(_.newDeclaration)

      val newWorkflowDefinitionElement: WorkflowDefinitionElement = {
        val withNewInputElements = a.definitionElement.copy(inputsSection = newInputsSection, graphElements = newGraphNodeSet)
        val identifierRenames = inputElementsWithUpstreams.map(ie => ie.original.name -> ie.newDeclaration.name).toMap
        withNewInputElements.renameIdentifiers(identifierRenames)
      }

      a.copy(definitionElement = newWorkflowDefinitionElement)

    } else a
  }
}
