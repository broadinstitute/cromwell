package wdl.transforms.base.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.KvPair

object AstNodeToKvPair {
  def astNodeToKvPair(implicit astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]): CheckedAtoB[GenericAstNode, KvPair] = CheckedAtoB.fromErrorOr {
    case a: GenericAst if a.getName == "ObjectKV" || a.getName == "MapLiteralKv" || a.getName == "RuntimeAttribute" =>
      val keyValidation: ErrorOr[String] = a.getAttributeAs[String]("key").toValidated


      val valueValidation: ErrorOr[ExpressionElement] = {
        val validation = a.getAttributeAs[ExpressionElement]("value").toValidated
        val forKeyContext: Option[String] = keyValidation.map(k => s"read value for key '$k'").toOption

        forKeyContext match {
          case Some(context) => validation.contextualizeErrors(context)
          case None => validation
        }

      }

      (keyValidation, valueValidation) mapN { (key, value) => KvPair(key, value) }
  }
}
