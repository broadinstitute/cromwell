package wdl.draft3.transforms.linking

import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.wdlom2wom.expression.WdlomWomExpression
import wdl.model.draft3.graph.UnlinkedValueConsumer.ops._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.WomExpressionMaker
import wdl.model.draft3.graph._
import wom.expression.WomExpression

package object expression {

  implicit val expressionElementSetUnlinkedValueConsumer: UnlinkedValueConsumer[Set[ExpressionElement]] = new UnlinkedValueConsumer[Set[ExpressionElement]] {
    override def consumedValueHooks(elements: Set[ExpressionElement]): Set[UnlinkedConsumedValueHook] = elements.flatMap { e: ExpressionElement => e.consumedValueHooks }
  }

  implicit val identifierLookupUnlinkedValueConsumer: UnlinkedValueConsumer[IdentifierLookup] = new UnlinkedValueConsumer[IdentifierLookup] {
    override def consumedValueHooks(a: IdentifierLookup): Set[UnlinkedConsumedValueHook] = Set(UnlinkedIdentifierHook(a.identifier))
  }

  implicit val objectLiteralUnlinkedValueConsumer: UnlinkedValueConsumer[ObjectLiteral] = new UnlinkedValueConsumer[ObjectLiteral] {
    override def consumedValueHooks(o: ObjectLiteral): Set[UnlinkedConsumedValueHook] = o.elements.values.toSet[ExpressionElement].consumedValueHooks
  }

  implicit val mapLiteralUnlinkedValueConsumer: UnlinkedValueConsumer[MapLiteral] = new UnlinkedValueConsumer[MapLiteral] {
    override def consumedValueHooks(m: MapLiteral): Set[UnlinkedConsumedValueHook] = m.elements.keys.toSet[ExpressionElement].consumedValueHooks ++ m.elements.values.toSet[ExpressionElement].consumedValueHooks
  }

  implicit val pairLiteralUnlinkedValueConsumer: UnlinkedValueConsumer[PairLiteral] = new UnlinkedValueConsumer[PairLiteral] {
    override def consumedValueHooks(p: PairLiteral): Set[UnlinkedConsumedValueHook] = p.left.consumedValueHooks ++ p.right.consumedValueHooks
  }

  implicit val arrayLiteralUnlinkedValueConsumer: UnlinkedValueConsumer[ArrayLiteral] = new UnlinkedValueConsumer[ArrayLiteral] {
    override def consumedValueHooks(a: ArrayLiteral): Set[UnlinkedConsumedValueHook] = a.elements.toSet.flatMap { e: ExpressionElement => e.consumedValueHooks }
  }

  implicit val identifierMemberAccessUnlinkedValueConsumer: UnlinkedValueConsumer[IdentifierMemberAccess] = new UnlinkedValueConsumer[IdentifierMemberAccess] {
    override def consumedValueHooks(a: IdentifierMemberAccess): Set[UnlinkedConsumedValueHook] = Set(UnlinkedCallOutputOrIdentifierAndMemberAccessHook(a.first, a.second))
  }

  implicit val expressionElementUnlinkedValueConsumer: UnlinkedValueConsumer[ExpressionElement] = new UnlinkedValueConsumer[ExpressionElement] {
    override def consumedValueHooks(a: ExpressionElement): Set[UnlinkedConsumedValueHook] = a match {
      case _: PrimitiveLiteralExpressionElement | _: StringLiteral => Set.empty
      case a: ObjectLiteral => a.consumedValueHooks
      case a: PairLiteral => a.consumedValueHooks
      case a: ArrayLiteral => a.consumedValueHooks
      case a: MapLiteral => a.consumedValueHooks

      case a: IdentifierLookup => a.consumedValueHooks
      case a: IdentifierMemberAccess => a.consumedValueHooks


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
