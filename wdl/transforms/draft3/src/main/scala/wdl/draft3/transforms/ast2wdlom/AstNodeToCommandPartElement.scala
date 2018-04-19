package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.validated._
import cats.syntax.either._
import cats.instances.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.draft3.parser.WdlParser._
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.{CommandPartElement, ExpressionElement}
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.elements.ExpressionElement.{PrimitiveLiteralExpressionElement, StringExpression, StringLiteral}

object AstNodeToCommandPartElement {
  def convert(a: AstNode): ErrorOr[CommandPartElement] = a match {
    case t: Terminal => astNodeToString(t).toValidated map StringCommandPartElement
    case a: Ast =>
      val expressionElementV: ErrorOr[ExpressionElement] = a.getAttributeAs[ExpressionElement]("expr").toValidated

      implicit val attributeConverter: CheckedAtoB[AstNode, (String, String)] = astNodeToAst andThen CheckedAtoB.fromErrorOr(convertAttributeKvp)

      val attributesV: ErrorOr[Map[String, String]] = a.getAttributeAsVector[(String, String)]("attributes").map(_.toMap).toValidated

      (expressionElementV, attributesV) mapN { (expressionElement, attributes) => PlaceholderCommandPartElement(expressionElement, attributes) }
    case other => s"Conversion for $other not supported".invalidNel
  }

  def convertAttributeKvp(a: Ast): ErrorOr[(String, String)] = {
    (a.getAttributeAs[String]("key").toValidated : ErrorOr[String],
      a.getAttributeAs[ExpressionElement]("value").toValidated : ErrorOr[ExpressionElement]) flatMapN { (key, value) =>
      val valueString: ErrorOr[String] = value match {
        case StringLiteral(literalValue) => literalValue.validNel
        case StringExpression(pieces) if pieces.length == 1 => pieces.head match {
          case StringLiteral(literalValue) => literalValue.validNel
          case other => s"Cannot use $other as a placeholder attribute. Must be a primitive literal".invalidNel
        }
        case PrimitiveLiteralExpressionElement(primitive) => primitive.valueString.validNel
        case other => s"Cannot use $other as a placeholder attribute. Must be a primitive literal".invalidNel
      }
      valueString.map(key -> _)
    }
  }
}
