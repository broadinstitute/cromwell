package wdl.transforms.base.wdlom2wom.graph

import cats.syntax.apply._
import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr.{ErrorOr, _}
import shapeless.Coproduct
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.model.draft3.graph.GraphElementValueConsumer.ops._
import wdl.model.draft3.graph.UnlinkedValueGenerator.ops._
import wdl.transforms.base.linking.graph._
import wdl.transforms.base.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition
import wdl.transforms.base.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.GraphLikeConvertInputs
import wdl.transforms.base.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.elements._
import wdl.model.draft3.graph._
import wdl.model.draft3.graph.expression.{FileEvaluator, TypeEvaluator, ValueEvaluator}
import wdl.shared.transforms.wdlom2wom.WomGraphMakerTools
import wom.callable.Callable.{InputDefinition, OverridableInputDefinitionWithDefault, OptionalInputDefinition, RequiredInputDefinition}
import wom.callable.{Callable, WorkflowDefinition}
import wom.graph.CallNode.{CallNodeBuilder, InputDefinitionFold, InputDefinitionPointer}
import wom.graph.GraphNode.GraphNodeSetter
import wom.graph.GraphNodePort.{ConnectedInputPort, InputPort, OutputPort}
import wom.graph._
import wom.graph.expression.{AnonymousExpressionNode, PlainAnonymousExpressionNode}
import wom.types.{WomAnyType, WomArrayType, WomType}

object ScatterElementToGraphNode {
  def convert(a: ScatterNodeMakerInputs)
             (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
              fileEvaluator: FileEvaluator[ExpressionElement],
              typeEvaluator: TypeEvaluator[ExpressionElement],
              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[GraphNode]] =
    if (a.convertNestedScatterToSubworkflow) {
      // Create a sub-workflow from the inner scatter.
      if (a.insideAnotherScatter) {
        convertInnerScatter(a)
      } else {
        convertOuterScatter(a)
      }
    } else {
      // do not do anything special with inner scatters.
      convertOuterScatter(a)
    }

  def convertOuterScatter(a: ScatterNodeMakerInputs)
                         (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
                          fileEvaluator: FileEvaluator[ExpressionElement],
                          typeEvaluator: TypeEvaluator[ExpressionElement],
                          valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[GraphNode]] = {
    val scatterExpression = a.node.scatterExpression
    val scatterVariableName = a.node.scatterVariableName
    val graphElements = a.node.graphElements

    val scatterWomExpressionV: ErrorOr[WdlomWomExpression] = WdlomWomExpression.make(scatterExpression, a.linkableValues)
    val scatterExpressionNodeValidation: ErrorOr[AnonymousExpressionNode] = scatterWomExpressionV flatMap { scatterWomExpression =>
      AnonymousExpressionNode.fromInputMapping(WomIdentifier(scatterVariableName), scatterWomExpression, a.linkablePorts, PlainAnonymousExpressionNode.apply)
    }

    val scatterVariableTypeValidation: ErrorOr[WomType] = scatterExpression.evaluateType(a.linkableValues) flatMap {
      case a: WomArrayType => a.memberType.validNel
      case WomAnyType => WomAnyType.validNel
      case other => s"Invalid type for scatter variable '$scatterVariableName': ${other.stableName}".invalidNel
    }

    final case class RequiredOuterPorts(valueGeneratorPorts: Map[String, OutputPort], completionPorts: Map[String, CallNode])

    val foundOuterGeneratorsValidation: ErrorOr[RequiredOuterPorts] = {
      val required: ErrorOr[Set[UnlinkedConsumedValueHook]] = graphElements.toList.traverse { element => element.graphElementConsumedValueHooks(a.availableTypeAliases, a.callables) }.map(_.toSet.flatten)
      val generated: ErrorOr[Set[GeneratedValueHandle]] = graphElements.toList.traverse { element => element.generatedValueHandles(a.availableTypeAliases, a.callables) }.map(_.toSet.flatten)

      def makeLink(hook: UnlinkedConsumedValueHook): (String, OutputPort) = {
        val name = a.linkableValues(hook).linkableName
        val port = a.linkablePorts(name)
        name -> port
      }

      (required, generated) mapN { (r, g) =>
        val requiredOuterValues = r collect {
          case hook@UnlinkedIdentifierHook(id) if id != scatterVariableName && !g.exists(_.linkableName == id) =>
            makeLink(hook)
          case hook@UnlinkedCallOutputOrIdentifierAndMemberAccessHook(first, second) if first != scatterVariableName && !g.exists(_.linkableName == first) && !g.exists(_.linkableName == s"$first.$second") =>
            makeLink(hook)
        }
        val requiredCompletionPorts = r collect {
          case UnlinkedAfterCallHook(upstreamCallName) if !g.exists {
            case GeneratedCallFinishedHandle(`upstreamCallName`) => true
            case _ => false
          } => upstreamCallName -> a.upstreamCalls.values.find(_.localName == upstreamCallName).get
        }

        RequiredOuterPorts(requiredOuterValues.toMap, requiredCompletionPorts.toMap)
      }
    }

    (scatterExpressionNodeValidation, scatterVariableTypeValidation, foundOuterGeneratorsValidation) flatMapN { (expressionNode, scatterVariableType, foundOuterGenerators) =>
      val womInnerGraphScatterVariableInput = ScatterVariableNode(WomIdentifier(scatterVariableName), expressionNode, scatterVariableType)
      val ogins: Set[GraphNode] = (foundOuterGenerators.valueGeneratorPorts.toList map { case (name: String, port: OutputPort) =>
        OuterGraphInputNode(WomIdentifier(name), port, preserveScatterIndex = false)
      }).toSet

      val graphLikeConvertInputs = GraphLikeConvertInputs(graphElements.toSet, ogins ++ Set(womInnerGraphScatterVariableInput), foundOuterGenerators.completionPorts, a.availableTypeAliases, a.workflowName,
                                                          insideAScatter = true,
                                                          convertNestedScatterToSubworkflow = a.convertNestedScatterToSubworkflow,
                                                          a.callables)
      val innerGraph: ErrorOr[Graph] = WorkflowDefinitionElementToWomWorkflowDefinition.convertGraphElements(graphLikeConvertInputs)

      innerGraph map { ig =>
        val withOutputs = WomGraphMakerTools.addDefaultOutputs(ig)
        val generatedAndNew = ScatterNode.scatterOverGraph(withOutputs, womInnerGraphScatterVariableInput)
        generatedAndNew.nodes
      }
    }
  }

  def convertInnerScatter(a: ScatterNodeMakerInputs)
                         (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
                          fileEvaluator: FileEvaluator[ExpressionElement],
                          typeEvaluator: TypeEvaluator[ExpressionElement],
                          valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[GraphNode]] = {

    val requiredOuterValuesValidation: ErrorOr[Set[UnlinkedConsumedValueHook]] = a.node.graphElementConsumedValueHooks(a.availableTypeAliases, a.callables)

    val subWorkflowInputsValidation: ErrorOr[Set[GraphNode]] = requiredOuterValuesValidation map { requiredOuterValues =>
      val requiredLinkableNames: Set[(String, WomType)] = requiredOuterValues map { hook => (a.linkableValues(hook).linkableName, a.linkableValues(hook).womType) }
      requiredLinkableNames map { case (name, womType) => RequiredGraphInputNode(WomIdentifier(name), womType, name) }
    }

    val subWorkflowGraphValidation: ErrorOr[Graph] = subWorkflowInputsValidation flatMap { subWorkflowInputs =>
      val graphLikeConvertInputs = GraphLikeConvertInputs(Set(a.node), subWorkflowInputs, Map.empty, a.availableTypeAliases, a.workflowName,
                                                          insideAScatter = false,
                                                          convertNestedScatterToSubworkflow = a.convertNestedScatterToSubworkflow,
                                                          a.callables)
      val subWorkflowGraph = WorkflowDefinitionElementToWomWorkflowDefinition.convertGraphElements(graphLikeConvertInputs)
      subWorkflowGraph map { WomGraphMakerTools.addDefaultOutputs(_) }
    }

    val subWorkflowDefinitionValidation = subWorkflowGraphValidation map { subWorkflowGraph =>
      WorkflowDefinition(a.node.scatterName,
                         subWorkflowGraph,
                         Map.empty,
                         Map.empty,
                         None // no lexical information, because this doesn't map to a source workflow
      )
    }

    val scatterableGraphValidation = subWorkflowDefinitionValidation map { subWorkflowDefinition =>
      val callNodeBuilder = new CallNodeBuilder()
      val graphNodeSetter = new GraphNodeSetter[CallNode]

      val unsatisfiedInputs = subWorkflowDefinition.inputs filter { i => !a.linkablePorts.contains(i.name) }
      val newInputNodes: Map[String, ExternalGraphInputNode] = (unsatisfiedInputs collect {
        case i: RequiredInputDefinition => i.name -> RequiredGraphInputNode(WomIdentifier(i.name), i.womType, i.name, Callable.InputDefinition.IdentityValueMapper)
        case i: OptionalInputDefinition => i.name -> OptionalGraphInputNode(WomIdentifier(i.name), i.womType, i.name, Callable.InputDefinition.IdentityValueMapper)
        case i: OverridableInputDefinitionWithDefault => i.name -> OptionalGraphInputNodeWithDefault(WomIdentifier(i.name), i.womType, i.default, i.name, Callable.InputDefinition.IdentityValueMapper)
      }).toMap

      val mappingAndPorts: List[((InputDefinition, InputDefinitionPointer), InputPort)] = subWorkflowDefinition.inputs map { i =>
        val port: OutputPort = a.linkablePorts.getOrElse(i.name, newInputNodes(i.name).singleOutputPort)
        val pointer = Coproduct[InputDefinitionPointer](port)
        (i -> pointer, ConnectedInputPort(i.name, i.womType, port, graphNodeSetter.get))
      }
      val mapping = mappingAndPorts.map(_._1)
      val inputPorts = mappingAndPorts.map(_._2).toSet
      val result = callNodeBuilder.build(WomIdentifier(a.node.scatterName),
                                         subWorkflowDefinition,
                                         InputDefinitionFold(mappings = mapping, callInputPorts = inputPorts),
                                         Set.empty,
                                         a.node.sourceLocation,
                                         (_, localName) => WomIdentifier(localName))
      graphNodeSetter._graphNode = result.node
      result.copy(newInputs = result.newInputs ++ newInputNodes.values)
    }

    scatterableGraphValidation map { scatterableGraph =>
      scatterableGraph.nodes
    }
  }
}

final case class ScatterNodeMakerInputs(node: ScatterElement,
                                        upstreamCalls: Map[String, CallNode],
                                        linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                        linkablePorts: Map[String, OutputPort],
                                        availableTypeAliases: Map[String, WomType],
                                        workflowName: String,
                                        insideAnotherScatter: Boolean,
                                        convertNestedScatterToSubworkflow: Boolean,
                                        callables: Map[String, Callable])
