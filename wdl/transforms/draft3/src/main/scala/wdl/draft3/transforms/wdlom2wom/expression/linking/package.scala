package wdl.draft3.transforms.wdlom2wom.expression

import wdl.draft3.transforms.wdlom2wom.linking._
import wdl.draft3.transforms.wdlom2wom.linking.UnlinkedValueConsumer.ops._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{IdentifierLookup, PrimitiveLiteralExpressionElement}

package object linking {

  implicit object expressionElementUnlinkedValueConsumer extends UnlinkedValueConsumer[ExpressionElement] {
    override def consumedValueNames(a: ExpressionElement): Set[UnlinkedConsumedValueName] = a match {
      case _: PrimitiveLiteralExpressionElement => Set.empty
      case id: IdentifierLookup => id.consumedValueNames
      // TODO fill in other expression types
      case other => throw new Exception(s"Cannot generate consumed values for ExpressionElement ${other.getClass.getSimpleName}")
    }
  }

  implicit object identifierLookupUnlinkedValueConsumer extends UnlinkedValueConsumer[IdentifierLookup] {
    override def consumedValueNames(a: IdentifierLookup): Set[UnlinkedConsumedValueName] = Set(UnlinkedIdentifierName(a.identifier))
  }
}
