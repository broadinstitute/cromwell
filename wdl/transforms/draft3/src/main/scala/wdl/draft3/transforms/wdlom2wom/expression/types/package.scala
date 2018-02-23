package wdl.draft3.transforms.wdlom2wom.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement.{IdentifierLookup, PrimitiveLiteralExpressionElement}
import wom.types.WomType

package object types {
  implicit object primitiveTypeEvaluator extends TypeEvaluator[PrimitiveLiteralExpressionElement] {
    override def evaluateType(a: PrimitiveLiteralExpressionElement, linkedValues: Set[LinkedConsumedValue]): ErrorOr[WomType] = a.value.womType.validNel
  }

  implicit object identifierLookupEvaluator extends TypeEvaluator[IdentifierLookup] {
    override def evaluateType(a: IdentifierLookup, linkedValues: Set[LinkedConsumedValue]): ErrorOr[WomType] = {
      linkedValues.collectFirst {
        case lcv if lcv.generatedValueName.linkableName == a.identifier => lcv.generatedType
      } match {
        case Some(womType) => womType.validNel
        case None => s"Type evaluation failure. No suitable input for identifier lookup '${a.identifier}' amongst {${linkedValues.map(_.generatedValueName.linkableName).mkString(", ")}}".invalidNel
      }
    }
  }

}
