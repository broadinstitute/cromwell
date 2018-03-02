package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.{Ast, AstNode}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wdl.model.draft3.elements.MetaValueElement
import wdl.model.draft3.elements.MetaValueElement.MetaKvPair

object AstNodeToMetaKvPair {
  def convert(astNode: AstNode): ErrorOr[MetaKvPair] = astNode match {
    case a: Ast if a.getName == "MetaKvPair" =>
      val keyValidation: ErrorOr[String] = a.getAttributeAs[String]("key").toValidated
      val valueValidation: ErrorOr[MetaValueElement] = a.getAttributeAs[MetaValueElement]("value").toValidated

      (keyValidation, valueValidation) mapN { (key, value) => MetaKvPair(key, value) }
  }
}
