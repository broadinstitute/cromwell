package wdl.draft3.transforms.linking.expression

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement.{IdentifierLookup, ObjectLiteral, PrimitiveLiteralExpressionElement}
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook, UnlinkedIdentifierHook}
import wdl.model.draft3.graph.expression.TypeEvaluator
import wom.types.{WomObjectType, WomType}

package object types {
  implicit val primitiveTypeEvaluator: TypeEvaluator[PrimitiveLiteralExpressionElement] = new TypeEvaluator[PrimitiveLiteralExpressionElement] {
    override def evaluateType(a: PrimitiveLiteralExpressionElement, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = a.value.womType.validNel
  }

  implicit val identifierLookupTypeEvaluator: TypeEvaluator[IdentifierLookup] = new TypeEvaluator[IdentifierLookup] {
    override def evaluateType(a: IdentifierLookup, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      linkedValues.collectFirst {
        case (UnlinkedIdentifierHook(id), gen) if a.identifier == id => gen.womType
      } match {
        case Some(womType) => womType.validNel
        case None => s"Type evaluation failure. No suitable input for identifier lookup '${a.identifier}' amongst {${linkedValues.map(_._2.linkableName).mkString(", ")}}".invalidNel
      }
    }
  }

  implicit val objectLiteralTypeEvaluator: TypeEvaluator[ObjectLiteral] = new TypeEvaluator[ObjectLiteral] {
    override def evaluateType(a: ObjectLiteral, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = WomObjectType.validNel
  }
}
