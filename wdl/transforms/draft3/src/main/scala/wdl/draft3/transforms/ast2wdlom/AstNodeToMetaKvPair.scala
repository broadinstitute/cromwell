package wdl.draft3.transforms.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.parser.WdlParser.{Ast, AstNode}
import wdl.draft3.transforms.ast2wdlom.EnhancedDraft3Ast._
import wom.callable.{MetaKvPair, MetaValueElement}

object AstNodeToMetaKvPair {
  def convert(astNode: AstNode): ErrorOr[MetaKvPair] = astNode match {
    case a: Ast if a.getName == "MetaKvPair" =>
      val keyValidation: ErrorOr[String] = a.getAttributeAs[String]("key").toValidated
      val valueValidation: ErrorOr[MetaValueElement] = a.getAttributeAs[MetaValueElement]("value").toValidated

      (keyValidation, valueValidation) mapN { (key, value) => MetaKvPair(key, value) }
  }
}
