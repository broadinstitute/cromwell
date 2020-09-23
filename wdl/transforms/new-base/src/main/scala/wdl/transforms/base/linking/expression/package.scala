package wdl.transforms.base.linking

import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import common.validation.ErrorOr._
import wdl.transforms.base.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.expression.{FileEvaluator, TypeEvaluator, ValueEvaluator, WomExpressionMaker}
import wdl.model.draft3.graph._
import wom.expression.WomExpression
import wom.types.WomType
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._

package object expression {

  implicit def expressionElementToWomExpression(implicit expressionValueConsumer: ExpressionValueConsumer[ExpressionElement],
                                                fileEvaluator: FileEvaluator[ExpressionElement],
                                                typeEvaluator: TypeEvaluator[ExpressionElement],
                                                valueEvaluator: ValueEvaluator[ExpressionElement]): WomExpressionMaker[ExpressionElement] = new WomExpressionMaker[ExpressionElement] {
    override def makeWomExpression(a: ExpressionElement,
                                   typeAliases: Map[String, WomType],
                                   consumedValueLookup: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomExpression] = {
      val consumedValueHooks = a.expressionConsumedValueHooks

      val neededLinkedValues = consumedValueHooks.toList.traverse {
        case c if consumedValueLookup.contains(c) => (c -> consumedValueLookup(c)).validNel
        case missing => s"Could not create WOM expression for '$a': Found no generated value for consumed value $missing".invalidNel
      }

      neededLinkedValues flatMap { lookup => WdlomWomExpression.make(a, lookup.toMap) }

    }
  }
}
