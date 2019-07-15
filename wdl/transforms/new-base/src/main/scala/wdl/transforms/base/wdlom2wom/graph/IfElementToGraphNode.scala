package wdl.transforms.base.wdlom2wom.graph

import cats.syntax.apply._
import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr.{ErrorOr, _}
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
import wom.callable.Callable
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.{AnonymousExpressionNode, PlainAnonymousExpressionNode}
import wom.types.{WomAnyType, WomBooleanType, WomType}

object IfElementToGraphNode {
  def convert(a: ConditionalNodeMakerInputs)
             (implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
              fileEvaluator: FileEvaluator[ExpressionElement],
              typeEvaluator: TypeEvaluator[ExpressionElement],
              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[GraphNode]] = {
    val conditionExpression = a.node.conditionExpression
    val graphElements = a.node.graphElements

    val conditionWomExpressionV: ErrorOr[WdlomWomExpression] = WdlomWomExpression.make(conditionExpression, a.linkableValues)
    val conditionExpressionNodeValidation: ErrorOr[AnonymousExpressionNode] = conditionWomExpressionV flatMap { conditionWomExpression =>
      AnonymousExpressionNode.fromInputMapping(WomIdentifier("if_condition"), conditionWomExpression, a.linkablePorts, PlainAnonymousExpressionNode.apply) }

    val conditionVariableTypeValidation: ErrorOr[Unit] = conditionExpression.evaluateType(a.linkableValues) flatMap {
      case WomBooleanType | WomAnyType => ().validNel
      case other => s"Invalid type for condition variable: ${other.stableName}".invalidNel
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
          case hook@UnlinkedIdentifierHook(id) if !g.exists(_.linkableName == id) => makeLink(hook)
          case hook@UnlinkedCallOutputOrIdentifierAndMemberAccessHook(first, second) if !g.exists(_.linkableName == first) && !g.exists(_.linkableName == s"$first.$second") => makeLink(hook)
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

    (conditionExpressionNodeValidation, conditionVariableTypeValidation, foundOuterGeneratorsValidation) flatMapN { (expressionNode, _, foundOuterGenerators) =>
      val ogins: Set[GraphNode] = (foundOuterGenerators.valueGeneratorPorts.toList map { case (name: String, port: OutputPort) =>
        OuterGraphInputNode(WomIdentifier(name), port, preserveScatterIndex = true)
      }).toSet

      val graphLikeConvertInputs = GraphLikeConvertInputs(graphElements.toSet, ogins, foundOuterGenerators.completionPorts, a.availableTypeAliases, a.workflowName,
                                                          insideAScatter = a.insideAnotherScatter,
                                                          convertNestedScatterToSubworkflow = a.convertNestedScatterToSubworkflow,
                                                          a.callables)
      val innerGraph: ErrorOr[Graph] = WorkflowDefinitionElementToWomWorkflowDefinition.convertGraphElements(graphLikeConvertInputs)

      innerGraph map { ig =>
        val withOutputs = WomGraphMakerTools.addDefaultOutputs(ig)
        val generatedAndNew = ConditionalNode.wireInConditional(withOutputs, expressionNode)
        generatedAndNew.nodes
      }
    }
  }
}

final case class ConditionalNodeMakerInputs(node: IfElement,
                                            upstreamCalls: Map[String, CallNode],
                                            linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                            linkablePorts: Map[String, OutputPort],
                                            availableTypeAliases: Map[String, WomType],
                                            workflowName: String,
                                            insideAnotherScatter: Boolean,
                                            convertNestedScatterToSubworkflow : Boolean,
                                            callables: Map[String, Callable])
