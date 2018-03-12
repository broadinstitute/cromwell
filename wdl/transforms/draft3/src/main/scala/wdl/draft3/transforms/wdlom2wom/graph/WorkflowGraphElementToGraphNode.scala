package wdl.draft3.transforms.wdlom2wom.graph

import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, _}
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.model.draft3.graph.expression.WomExpressionMaker.ops._
import wdl.draft3.transforms.linking.typemakers._
import wdl.draft3.transforms.linking.expression._
import wdl.draft3.transforms.linking.expression.types._
import wdl.draft3.transforms.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition
import wdl.draft3.transforms.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.GraphLikeConvertInputs
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.elements._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.shared.transforms.wdlom2wom.WomGraphMakerTools
import wom.expression.WomExpression
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.{AnonymousExpressionNode, ExposedExpressionNode, PlainAnonymousExpressionNode}
import wom.types.{WomArrayType, WomType}

object WorkflowGraphElementToGraphNode {
  def convert(a: GraphNodeMakerInputs): ErrorOr[Set[GraphNode]] = a.node match {
    case InputDeclarationElement(typeElement, name, None) =>
      typeElement.determineWomType(a.availableTypeAliases) map { womType =>
        Set(RequiredGraphInputNode(WomIdentifier(name), womType, s"${a.workflowName}.$name"))
      }
    case DeclarationElement(typeElement, name, Some(expr)) =>
      val womExprValidation: ErrorOr[WomExpression] = expr.makeWomExpression(a.availableTypeAliases, a.linkableValues)
      val womTypeValidation: ErrorOr[WomType] = typeElement.determineWomType(a.availableTypeAliases)

      (womExprValidation, womTypeValidation) flatMapN { (womExpr, womType) =>
        a.node match {
          case _: InputDeclarationElement =>
            Set[GraphNode](OptionalGraphInputNodeWithDefault.apply(WomIdentifier(name), womType, womExpr, name)).validNel
          case _: IntermediateValueDeclarationElement =>
            ExposedExpressionNode.fromInputMapping(WomIdentifier(name), womExpr, womType, a.linkablePorts) map { Set(_) }
          case _: OutputDeclarationElement =>
            ExpressionBasedGraphOutputNode.fromInputMapping(WomIdentifier(s"${a.workflowName}.$name"), womExpr, womType, a.linkablePorts) map {Set(_)}
        }
      }

    case ScatterElement(scatterExpression, scatterVariableName, graphElements) =>
      val scatterWomExpression: WdlomWomExpression = WdlomWomExpression(scatterExpression, a.linkableValues)
      val scatterExpressionNodeValidation: ErrorOr[AnonymousExpressionNode] = AnonymousExpressionNode.fromInputMapping(WomIdentifier(scatterVariableName), scatterWomExpression, a.linkablePorts, PlainAnonymousExpressionNode.apply)

      val scatterVariableTypeValidation: ErrorOr[WomType] = scatterExpression.evaluateType(a.linkableValues) flatMap {
        case a: WomArrayType => a.memberType.validNel
        case _ => "Invalid type for scatter variable: ".invalidNel
      }

      (scatterExpressionNodeValidation, scatterVariableTypeValidation) flatMapN { (expressionNode, scatterVariableType) =>
        val womInnerGraphScatterVariableInput = ScatterVariableNode(WomIdentifier(scatterVariableName), expressionNode, scatterVariableType)

        val graphLikeConvertInputs = GraphLikeConvertInputs(graphElements.toSet, Set(womInnerGraphScatterVariableInput), a.availableTypeAliases, a.workflowName)
        val innerGraph: ErrorOr[Graph] = WorkflowDefinitionElementToWomWorkflowDefinition.convertGraphElements(graphLikeConvertInputs)

        innerGraph map { ig =>
          val withOutputs = WomGraphMakerTools.addDefaultOutputs(ig)
          val generatedAndNew = ScatterNode.scatterOverGraph(withOutputs, womInnerGraphScatterVariableInput)
          generatedAndNew.nodes
        }
      }
  }
}

final case class GraphNodeMakerInputs(node: WorkflowGraphElement,
                                      linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                      linkablePorts: Map[String, OutputPort],
                                      availableTypeAliases: Map[String, WomType],
                                      workflowName: String)
