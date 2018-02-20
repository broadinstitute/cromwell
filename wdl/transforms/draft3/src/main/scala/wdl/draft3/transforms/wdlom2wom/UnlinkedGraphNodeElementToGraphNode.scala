package wdl.draft3.transforms.wdlom2wom

import cats.syntax.either._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.wdlom2wom.TypeElementToType.TypeElementToWomTypeParameters
import wdl.model.draft3.elements.{DeclarationElement, InputDeclarationElement, IntermediateValueDeclarationElement, OutputDeclarationElement}
import wom.expression.WomExpression
import wom.graph.GraphNodePort.OutputPort
import wom.graph.expression.ExposedExpressionNode
import wom.graph._
import wom.types.WomType

object UnlinkedGraphNodeElementToGraphNode {
  def convert(a: GraphNodeMakerInputs): ErrorOr[GraphNode] = a.node.fromElement match {
    case InputDeclarationElement(typeElement, name, None) =>
      (TypeElementToType.convert(TypeElementToWomTypeParameters(typeElement)) map { womType =>
        RequiredGraphInputNode(WomIdentifier(name), womType, name)
      }).toValidated
    case DeclarationElement(typeElement, name, Some(expr)) =>
      val womExprValidation: ErrorOr[WomExpression] = ExpressionElementToWomExpression.convert(WomExpressionMakerInputs(expr, a.linkableValues))
      val womTypeValidation: ErrorOr[WomType] = TypeElementToType.convert(TypeElementToWomTypeParameters(typeElement)).toValidated

      (womExprValidation, womTypeValidation) flatMapN { (womExpr, womType) =>
        a.node.fromElement match {
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

final case class GraphNodeMakerInputs(node: UnlinkedGraphNode, linkableValues: Map[String, OutputPort], workflowName: String)
