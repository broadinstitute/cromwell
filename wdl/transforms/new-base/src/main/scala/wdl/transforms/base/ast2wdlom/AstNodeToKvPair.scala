package wdl.transforms.base.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.KvPair

object AstNodeToKvPair {
  def astNodeToKvPair(implicit
    astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]
  ): CheckedAtoB[GenericAstNode, KvPair] = CheckedAtoB.fromErrorOr {
    case a: GenericAst if a.getName == "ObjectKV" || a.getName == "MapLiteralKv" || a.getName == "RuntimeAttribute" =>
      validate(a.getAttributeAs[String]("key"), a.getAttributeAs[ExpressionElement]("value"))
  }

  def validate(key: Checked[String], value: Checked[ExpressionElement]): ErrorOr[KvPair] = {
    val keyValidation = key.toValidated
    val valueValidation = value.toValidated
    val forKeyContext: Option[String] = keyValidation.map(k => s"read value for key '$k'").toOption

    forKeyContext match {
      case Some(context) => valueValidation.contextualizeErrors(context)
      case None => valueValidation
    }
    (keyValidation, valueValidation) mapN { (key, value) => KvPair(key, value) }
  }
}
