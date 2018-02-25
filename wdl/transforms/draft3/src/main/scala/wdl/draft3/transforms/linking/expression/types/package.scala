package wdl.draft3.transforms.linking.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement.{IdentifierLookup, PrimitiveLiteralExpressionElement}
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook, UnlinkedIdentifierHook}
import wdl.model.draft3.graph.expression.TypeEvaluator
import wom.types.WomType

package object types {
  implicit object primitiveTypeEvaluator extends TypeEvaluator[PrimitiveLiteralExpressionElement] {
    override def evaluateType(a: PrimitiveLiteralExpressionElement, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = a.value.womType.validNel
  }

  implicit object identifierLookupTypeEvaluator extends TypeEvaluator[IdentifierLookup] {
    override def evaluateType(a: IdentifierLookup, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      linkedValues.collectFirst {
        case (UnlinkedIdentifierHook(id), gen) if a.identifier == id => gen.womType
      } match {
        case Some(womType) => womType.validNel
        case None => s"Type evaluation failure. No suitable input for identifier lookup '${a.identifier}' amongst {${linkedValues.map(_._2.linkableName).mkString(", ")}}".invalidNel
      }
    }
  }
}
