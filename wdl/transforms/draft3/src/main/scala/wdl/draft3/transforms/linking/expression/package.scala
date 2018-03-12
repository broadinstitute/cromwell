package wdl.draft3.transforms.linking

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.expression.WomExpressionMaker
import wdl.model.draft3.graph._
import wom.expression.WomExpression
import wom.types.WomType
import wdl.model.draft3.graph.ExpressionValueConsumer.ops._
import wdl.draft3.transforms.linking.expression.consumed._

package object expression {

  implicit val expressionElementToWomExpression: WomExpressionMaker[ExpressionElement] = new WomExpressionMaker[ExpressionElement] {
    override def makeWomExpression(a: ExpressionElement,
                                   typeAliases: Map[String, WomType],
                                   consumedValueLookup: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomExpression] = {
      val consumedValueHooks = a.expressionConsumedValueHooks

      val neededLinkedValues = consumedValueHooks.toList.traverse[ErrorOr, (UnlinkedConsumedValueHook, GeneratedValueHandle)] {
        case c if consumedValueLookup.contains(c) => (c -> consumedValueLookup(c)).validNel
        case missing => s"Could not create WOM expression for '$a': Found no generated value for consumed value $missing".invalidNel
      }

      neededLinkedValues map { lookup => WdlomWomExpression(a, lookup.toMap) }

    }
  }
}
