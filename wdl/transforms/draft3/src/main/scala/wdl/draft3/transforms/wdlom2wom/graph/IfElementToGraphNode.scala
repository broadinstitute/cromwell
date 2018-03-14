package wdl.draft3.transforms.wdlom2wom.graph

import cats.syntax.apply._
import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr.{ErrorOr, _}
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
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.{AnonymousExpressionNode, PlainAnonymousExpressionNode}
import wom.types.{WomBooleanType, WomType}

object IfElementToGraphNode {
  def convert(a: ConditionalNodeMakerInputs): ErrorOr[Set[GraphNode]] = {
    val conditionExpression = a.node.conditionExpression
    val graphElements = a.node.graphElements

    val conditionWomExpression: WdlomWomExpression = WdlomWomExpression(conditionExpression, a.linkableValues)
    val conditionExpressionNodeValidation: ErrorOr[AnonymousExpressionNode] = AnonymousExpressionNode.fromInputMapping(WomIdentifier("if_condition"), conditionWomExpression, a.linkablePorts, PlainAnonymousExpressionNode.apply)

    val scatterVariableTypeValidation: ErrorOr[Unit] = conditionExpression.evaluateType(a.linkableValues) flatMap {
      case WomBooleanType => ().validNel
      case other => s"Invalid type for condition variable: ${other.toDisplayString}".invalidNel
    }

    val requiredOuterValuesValidation: ErrorOr[Set[String]] = {
      val required: ErrorOr[Set[UnlinkedConsumedValueHook]] = (graphElements.toList.traverse { element => element.graphElementConsumedValueHooks(a.availableTypeAliases) }).map(_.toSet.flatten)
      val generated: ErrorOr[Set[GeneratedValueHandle]] = (graphElements.toList.traverse { element => element.generatedValueHandles(a.availableTypeAliases) }).map(_.toSet.flatten)

      (required, generated) mapN { (r, g) => r collect {
        case hook @ UnlinkedIdentifierHook(id) if !g.exists(_.linkableName == id) => a.linkableValues(hook).linkableName
        case hook @ UnlinkedCallOutputOrIdentifierAndMemberAccessHook(first, second) if !g.exists(_.linkableName == first) && !g.exists(_.linkableName == s"$first.$second") => a.linkableValues(hook).linkableName
      }}
    }

    val foundOuterGeneratorsValidation: ErrorOr[Map[String, OutputPort]] = requiredOuterValuesValidation map { requiredOuterValues =>
      (requiredOuterValues map { name =>
        name -> a.linkablePorts(name)
      }).toMap
    }

    (conditionExpressionNodeValidation, scatterVariableTypeValidation, foundOuterGeneratorsValidation) flatMapN { (expressionNode, _, foundOuterGenerators) =>
      val ogins: Set[GraphNode] = (foundOuterGenerators.toList map { case (name: String, port: OutputPort) =>
        OuterGraphInputNode(WomIdentifier(name), port, preserveScatterIndex = false)
      }).toSet

      val graphLikeConvertInputs = GraphLikeConvertInputs(graphElements.toSet, ogins, a.availableTypeAliases, a.workflowName, insideAScatter = true)
      val innerGraph: ErrorOr[Graph] = WorkflowDefinitionElementToWomWorkflowDefinition.convertGraphElements(graphLikeConvertInputs)

      innerGraph map { ig =>
        val withOutputs = WomGraphMakerTools.addDefaultOutputs(ig, compoundCallIdentifiers = false)
        val generatedAndNew = ConditionalNode.wireInConditional(withOutputs, expressionNode)
        generatedAndNew.nodes
      }
    }
  }
}

final case class ConditionalNodeMakerInputs(node: IfElement,
                                            linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                            linkablePorts: Map[String, OutputPort],
                                            availableTypeAliases: Map[String, WomType],
                                            workflowName: String,
                                            insideAnotherScatter: Boolean)
