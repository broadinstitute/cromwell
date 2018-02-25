package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.{Ast, AstNode}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.KvPair

object AstNodeToKvPair {
  def convert(astNode: AstNode): ErrorOr[KvPair] = astNode match {
    case a: Ast if a.getName == "ObjectKV" || a.getName == "MapLiteralKv" =>
      val keyValidation: ErrorOr[String] = a.getAttributeAs[String]("key").toValidated
      val valueValidation: ErrorOr[ExpressionElement] = a.getAttributeAs[ExpressionElement]("value").toValidated

      (keyValidation, valueValidation) mapN { (key, value) => KvPair(key, value) }
  }
}
