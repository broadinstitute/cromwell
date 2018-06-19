package wdl.draft3.transforms.wdlom2wom.graph

import cats.syntax.apply._
import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr.{ErrorOr, _}
import shapeless.Coproduct
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.model.draft3.graph.GraphElementValueConsumer.ops._
import wdl.model.draft3.graph.UnlinkedValueGenerator.ops._
import wdl.draft3.transforms.linking.expression.types._
import wdl.draft3.transforms.linking.graph._
import wdl.draft3.transforms.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition
import wdl.draft3.transforms.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.GraphLikeConvertInputs
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.elements._
import wdl.model.draft3.graph._
import wdl.shared.transforms.wdlom2wom.WomGraphMakerTools
import wom.callable.Callable.{InputDefinition, InputDefinitionWithDefault, OptionalInputDefinition, RequiredInputDefinition}
import wom.callable.{Callable, WorkflowDefinition}
import wom.graph.CallNode.{CallNodeBuilder, InputDefinitionFold, InputDefinitionPointer}
import wom.graph.GraphNode.GraphNodeSetter
import wom.graph.GraphNodePort.{ConnectedInputPort, InputPort, OutputPort}
import wom.graph._
import wom.graph.expression.{AnonymousExpressionNode, PlainAnonymousExpressionNode}
import wom.types.{WomAnyType, WomArrayType, WomType}

object ScatterElementToGraphNode {
  def convert(a: ScatterNodeMakerInputs): ErrorOr[Set[GraphNode]] = if (a.insideAnotherScatter) {
    convertInnerScatter(a)
  } else {
    convertOuterScatter(a)
  }

  def convertOuterScatter(a: ScatterNodeMakerInputs): ErrorOr[Set[GraphNode]] = {
    val scatterExpression = a.node.scatterExpression
    val scatterVariableName = a.node.scatterVariableName
    val graphElements = a.node.graphElements

    val scatterWomExpression: WdlomWomExpression = WdlomWomExpression(scatterExpression, a.linkableValues)
    val scatterExpressionNodeValidation: ErrorOr[AnonymousExpressionNode] = AnonymousExpressionNode.fromInputMapping(WomIdentifier(scatterVariableName), scatterWomExpression, a.linkablePorts, PlainAnonymousExpressionNode.apply)

    val scatterVariableTypeValidation: ErrorOr[WomType] = scatterExpression.evaluateType(a.linkableValues) flatMap {
      case a: WomArrayType => a.memberType.validNel
      case WomAnyType => WomAnyType.validNel
      case other => s"Invalid type for scatter variable '$scatterVariableName': ${other.toDisplayString}".invalidNel
    }

    val requiredOuterValuesValidation: ErrorOr[Set[String]] = {
      val required: ErrorOr[Set[UnlinkedConsumedValueHook]] = graphElements.toList.traverse { element => element.graphElementConsumedValueHooks(a.availableTypeAliases, a.callables) }.map(_.toSet.flatten)
      val generated: ErrorOr[Set[GeneratedValueHandle]] = graphElements.toList.traverse { element => element.generatedValueHandles(a.availableTypeAliases, a.callables) }.map(_.toSet.flatten)

      (required, generated) mapN { (r, g) => r collect {
        case hook @ UnlinkedIdentifierHook(id) if id != scatterVariableName && !g.exists(_.linkableName == id) => a.linkableValues(hook).linkableName
        case hook @ UnlinkedCallOutputOrIdentifierAndMemberAccessHook(first, second)
          // NB: mirrors logic in scatterElementUnlinkedValueConsumer
          if !g.exists(_.linkableName == first) && !g.exists(_.linkableName == s"$first.$second") && (first != scatterVariableName) => a.linkableValues(hook).linkableName
      }}
    }

    val foundOuterGeneratorsValidation: ErrorOr[Map[String, OutputPort]] = requiredOuterValuesValidation map { requiredOuterValues =>
      (requiredOuterValues map { name =>
        name -> a.linkablePorts(name)
      }).toMap
    }

    (scatterExpressionNodeValidation, scatterVariableTypeValidation, foundOuterGeneratorsValidation) flatMapN { (expressionNode, scatterVariableType, foundOuterGenerators) =>
      val womInnerGraphScatterVariableInput = ScatterVariableNode(WomIdentifier(scatterVariableName), expressionNode, scatterVariableType)
      val ogins: Set[GraphNode] = (foundOuterGenerators.toList map { case (name: String, port: OutputPort) =>
        OuterGraphInputNode(WomIdentifier(name), port, preserveScatterIndex = false)
      }).toSet

      val graphLikeConvertInputs = GraphLikeConvertInputs(graphElements.toSet, ogins ++ Set(womInnerGraphScatterVariableInput), a.availableTypeAliases, a.workflowName, insideAScatter = true, a.callables)
      val innerGraph: ErrorOr[Graph] = WorkflowDefinitionElementToWomWorkflowDefinition.convertGraphElements(graphLikeConvertInputs)

      innerGraph map { ig =>
        val withOutputs = WomGraphMakerTools.addDefaultOutputs(ig)
        val generatedAndNew = ScatterNode.scatterOverGraph(withOutputs, womInnerGraphScatterVariableInput)
        generatedAndNew.nodes
      }
    }
  }

  def convertInnerScatter(a: ScatterNodeMakerInputs): ErrorOr[Set[GraphNode]] = {

    val requiredOuterValuesValidation: ErrorOr[Set[UnlinkedConsumedValueHook]] = a.node.graphElementConsumedValueHooks(a.availableTypeAliases, a.callables)

    val subWorkflowInputsValidation: ErrorOr[Set[GraphNode]] = requiredOuterValuesValidation map { requiredOuterValues =>
      val requiredLinkableNames: Set[(String, WomType)] = requiredOuterValues map { hook => (a.linkableValues(hook).linkableName, a.linkableValues(hook).womType) }
      requiredLinkableNames map { case (name, womType) => RequiredGraphInputNode(WomIdentifier(name), womType, name) }
    }

    val subWorkflowGraphValidation: ErrorOr[Graph] = subWorkflowInputsValidation flatMap { subWorkflowInputs =>
      val graphLikeConvertInputs = GraphLikeConvertInputs(Set(a.node), subWorkflowInputs, a.availableTypeAliases, a.workflowName, insideAScatter = false, a.callables)
      val subWorkflowGraph = WorkflowDefinitionElementToWomWorkflowDefinition.convertGraphElements(graphLikeConvertInputs)
      subWorkflowGraph map { WomGraphMakerTools.addDefaultOutputs(_) }
    }

    val subWorkflowDefinitionValidation = subWorkflowGraphValidation map { subWorkflowGraph => WorkflowDefinition(a.node.scatterName, subWorkflowGraph, Map.empty, Map.empty) }

    val scatterableGraphValidation = subWorkflowDefinitionValidation map { subWorkflowDefinition =>
      val callNodeBuilder = new CallNodeBuilder()
      val graphNodeSetter = new GraphNodeSetter[CallNode]

      val unsatisfiedInputs = subWorkflowDefinition.inputs filter { i => !a.linkablePorts.contains(i.name) }
      val newInputNodes: Map[String, ExternalGraphInputNode] = (unsatisfiedInputs collect {
        case i: RequiredInputDefinition => i.name -> RequiredGraphInputNode(WomIdentifier(i.name), i.womType, i.name, Callable.InputDefinition.IdentityValueMapper)
        case i: OptionalInputDefinition => i.name -> OptionalGraphInputNode(WomIdentifier(i.name), i.womType, i.name, Callable.InputDefinition.IdentityValueMapper)
        case i: InputDefinitionWithDefault => i.name -> OptionalGraphInputNodeWithDefault(WomIdentifier(i.name), i.womType, i.default, i.name, Callable.InputDefinition.IdentityValueMapper)
      }).toMap

      val mappingAndPorts: List[((InputDefinition, InputDefinitionPointer), InputPort)] = subWorkflowDefinition.inputs map { i =>
        val port: OutputPort = a.linkablePorts.getOrElse(i.name, newInputNodes(i.name).singleOutputPort)
        val pointer = Coproduct[InputDefinitionPointer](port)
        (i -> pointer, ConnectedInputPort(i.name, i.womType, port, graphNodeSetter.get))
      }
      val mapping = mappingAndPorts.map(_._1)
      val inputPorts = mappingAndPorts.map(_._2).toSet
      val result = callNodeBuilder.build(WomIdentifier(a.node.scatterName), subWorkflowDefinition, InputDefinitionFold(mappings = mapping, callInputPorts = inputPorts), (_, localName) => WomIdentifier(localName))
      graphNodeSetter._graphNode = result.node
      result.copy(newInputs = result.newInputs ++ newInputNodes.values)
    }

    scatterableGraphValidation map { scatterableGraph =>
      scatterableGraph.nodes
    }
  }
}

final case class ScatterNodeMakerInputs(node: ScatterElement,
                                        linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                        linkablePorts: Map[String, OutputPort],
                                        availableTypeAliases: Map[String, WomType],
                                        workflowName: String,
                                        insideAnotherScatter: Boolean,
                                        callables: Map[String, Callable])
