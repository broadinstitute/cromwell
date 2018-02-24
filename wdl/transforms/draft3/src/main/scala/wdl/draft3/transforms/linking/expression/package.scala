package wdl.draft3.transforms.linking

import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._

import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.graph.UnlinkedValueConsumer.ops._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{IdentifierLookup, PrimitiveLiteralExpressionElement}
import wdl.model.draft3.graph.expression.WomExpressionMaker
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook, UnlinkedIdentifierHook, UnlinkedValueConsumer}
import wom.expression.WomExpression

package object expression {

  implicit object identifierLookupUnlinkedValueConsumer extends UnlinkedValueConsumer[IdentifierLookup] {
    override def consumedValueHooks(a: IdentifierLookup): Set[UnlinkedConsumedValueHook] = Set(UnlinkedIdentifierHook(a.identifier))
  }

  implicit object expressionElementUnlinkedValueConsumer extends UnlinkedValueConsumer[ExpressionElement] {
    override def consumedValueHooks(a: ExpressionElement): Set[UnlinkedConsumedValueHook] = a match {
      case _: PrimitiveLiteralExpressionElement => Set.empty
      case id: IdentifierLookup => id.consumedValueHooks
      // TODO fill in other expression types
      case other => throw new Exception(s"Cannot generate consumed values for ExpressionElement ${other.getClass.getSimpleName}")
    }
  }

  implicit object expressionElementToWomExpression extends WomExpressionMaker[ExpressionElement] {
    override def makeWomExpression(a: ExpressionElement,
                                   consumedValueLookup: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomExpression] = {
      val neededLinkedValues = a.consumedValueHooks.toList.traverse[ErrorOr, (UnlinkedConsumedValueHook, GeneratedValueHandle)] {
        case c if consumedValueLookup.contains(c) => (c -> consumedValueLookup(c)).validNel
        case missing => s"Could not create WOM expression for '$a': Found no generated value for consumed value $missing".invalidNel
      }

      neededLinkedValues map { lookup => WdlomWomExpression(a, lookup.toMap) }
    }
  }
}
