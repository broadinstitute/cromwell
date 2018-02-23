package wdl.draft3.transforms.wdlom2wom.graph

import cats.syntax.either._
import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, _}
import wdl.draft3.transforms.wdlom2wom.TypeElementToWomType
import wdl.draft3.transforms.wdlom2wom.TypeElementToWomType.TypeElementToWomTypeParameters
import wdl.draft3.transforms.wdlom2wom.expression.{ExpressionElementToWomExpression, WomExpressionMakerInputs}
import wdl.model.draft3.elements._
import wom.expression.WomExpression
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.ExposedExpressionNode
import wom.types.WomType

object WorkflowGraphElementToGraphNode {
  def convert(a: GraphNodeMakerInputs): ErrorOr[GraphNode] = a.node match {
    case InputDeclarationElement(typeElement, name, None) =>
      (TypeElementToWomType.convert(TypeElementToWomTypeParameters(typeElement)) map { womType =>
        RequiredGraphInputNode(WomIdentifier(name), womType, name)
      }).toValidated
    case DeclarationElement(typeElement, name, Some(expr)) =>
      val womExprValidation: ErrorOr[WomExpression] = ExpressionElementToWomExpression.convert(WomExpressionMakerInputs(expr, a.linkableValues))
      val womTypeValidation: ErrorOr[WomType] = TypeElementToWomType.convert(TypeElementToWomTypeParameters(typeElement)).toValidated

      (womExprValidation, womTypeValidation) flatMapN { (womExpr, womType) =>
        a.node match {
          case _: InputDeclarationElement =>
            OptionalGraphInputNodeWithDefault.apply(WomIdentifier(name), womType, womExpr, name).validNel : ErrorOr[GraphNode]
          case _: IntermediateValueDeclarationElement =>
            ExposedExpressionNode.fromInputMapping(WomIdentifier(name), womExpr, womType, a.linkableValues)
          case _: OutputDeclarationElement =>
            ExpressionBasedGraphOutputNode.fromInputMapping(WomIdentifier(s"${a.workflowName}.$name"), womExpr, womType, a.linkableValues)
        }

      }
  }
}

final case class GraphNodeMakerInputs(node: WorkflowGraphElement, linkableValues: Map[String, OutputPort], workflowName: String)
