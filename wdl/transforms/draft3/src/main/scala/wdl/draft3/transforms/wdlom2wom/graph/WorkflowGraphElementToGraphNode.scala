package wdl.draft3.transforms.wdlom2wom.graph

import cats.syntax.validated._
import cats.syntax.traverse._
import cats.instances.list._
import common.validation.ErrorOr.{ErrorOr, _}
import wdl.model.draft3.graph.expression.WomTypeMaker.ops._
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wdl.model.draft3.graph.GraphElementValueConsumer.ops._
import wdl.model.draft3.graph.expression.WomExpressionMaker.ops._
import wdl.draft3.transforms.linking.typemakers._
import wdl.draft3.transforms.linking.expression._
import wdl.draft3.transforms.linking.expression.types._
import wdl.draft3.transforms.linking.graph._
import wdl.draft3.transforms.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition
import wdl.draft3.transforms.wdlom2wom.WorkflowDefinitionElementToWomWorkflowDefinition.GraphLikeConvertInputs
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.elements._
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedCallOutputOrIdentifierAndMemberAccessHook, UnlinkedConsumedValueHook, UnlinkedIdentifierHook}
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

    case se: ScatterElement =>
      val scatterMakerInputs = ScatterNodeMakerInputs(se, a.linkableValues, a.linkablePorts, a.availableTypeAliases, a.workflowName, a.insideAScatter)
      ScatterElementToGraphNode.convert(scatterMakerInputs)
  }
}

final case class GraphNodeMakerInputs(node: WorkflowGraphElement,
                                      linkableValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle],
                                      linkablePorts: Map[String, OutputPort],
                                      availableTypeAliases: Map[String, WomType],
                                      workflowName: String,
                                      insideAScatter: Boolean)
