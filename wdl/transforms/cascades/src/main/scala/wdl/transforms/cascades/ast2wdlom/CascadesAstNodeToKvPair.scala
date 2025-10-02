package wdl.transforms.cascades.ast2wdlom

import cats.syntax.apply._
import cats.syntax.either._
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr.{ErrorOr, _}
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.KvPair
import wdl.transforms.base.ast2wdlom.{GenericAst, GenericAstNode}

object CascadesAstNodeToKvPair {
  def astNodeToKvPair(implicit
    astNodeToExpressionElement: CheckedAtoB[GenericAstNode, ExpressionElement]
  ): CheckedAtoB[GenericAstNode, KvPair] = CheckedAtoB.fromErrorOr {
    case i: GenericAst if i.getName == "CallInputKV" =>
      // For CallInputKV, we accept either a key/value pair (`foo=foo`) or a plain identifier (`foo`).
      // In the latter case, "foo" is both the key and the name of an in-scope variable whose value should
      // be passed through. Due to limitations in the parser, both map to the same AST structure,
      // with the "value" attribute missing when a passthrough is intended.
      val valueKeyName = if (i.getAttributes.contains("value")) "value" else "key"
      validateKvPair(i.getAttributeAs[String]("key"), i.getAttributeAs[ExpressionElement](valueKeyName))
    case a: GenericAst if a.getName == "ObjectKV" || a.getName == "MapLiteralKv" || a.getName == "RuntimeAttribute" =>
      validateKvPair(a.getAttributeAs[String]("key"), a.getAttributeAs[ExpressionElement]("value"))
  }

  private def validateKvPair(key: Checked[String], value: Checked[ExpressionElement]): ErrorOr[KvPair] = {
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
