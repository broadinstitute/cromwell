package wdl.transforms.base.linking.expression.types

import cats.data.NonEmptyList
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import common.validation.Validation._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph._
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wom.types._

object LookupEvaluators {

  implicit val identifierLookupTypeEvaluator: TypeEvaluator[IdentifierLookup] = new TypeEvaluator[IdentifierLookup] {
    override def evaluateType(a: IdentifierLookup, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {
      linkedValues.collectFirst {
        case (UnlinkedIdentifierHook(id), gen) if a.identifier == id => gen.womType
      } match {
        case Some(womType) => womType.validNel
        case None => s"Type evaluation failure. No suitable type found for identifier lookup '${a.identifier}' amongst {${linkedValues.map(_._2.linkableName).mkString(", ")}}".invalidNel
      }
    }
  }

  implicit val expressionMemberAccessEvaluator: TypeEvaluator[ExpressionMemberAccess] = new TypeEvaluator[ExpressionMemberAccess] {
    override def evaluateType(a: ExpressionMemberAccess, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {
      val baseType = a.expression.evaluateType(linkedValues)
      baseType flatMap { doLookup(_, a.memberAccessTail) }
    }
  }

  implicit val identifierMemberAccessEvaluator: TypeEvaluator[IdentifierMemberAccess] = new TypeEvaluator[IdentifierMemberAccess] {
    override def evaluateType(a: IdentifierMemberAccess, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {
      val generatedValueHandle = linkedValues.get(UnlinkedCallOutputOrIdentifierAndMemberAccessHook(a.first, a.second))

      generatedValueHandle match {
        case Some(GeneratedIdentifierValueHandle(a.first, womType)) => doLookup(womType, NonEmptyList(a.second, a.memberAccessTail.toList))
        case Some(GeneratedCallOutputValueHandle(a.first, a.second, womType)) => NonEmptyList.fromList(a.memberAccessTail.toList) match {
          case Some(tailList) => doLookup(womType, tailList)
          case None => womType.validNel
        }
        case _ => s"Type evaluation failure. No suitable type found for identifier lookup '${a.first}' or '${a.first}.${a.second}' amongst {${linkedValues.map(_._2.linkableName).mkString(", ")}}".invalidNel
      }
    }
  }

  implicit val indexAccessTypeEvaluator: TypeEvaluator[IndexAccess] = new TypeEvaluator[IndexAccess] {
    override def evaluateType(a: IndexAccess, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle])
                             (implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] = {
      (a.expressionElement.evaluateType(linkedValues), a.index.evaluateType(linkedValues), a.index.validNel) flatMapN {
        case (a: WomArrayType, WomIntegerType, _) => a.memberType.validNel
        case (WomMapType(keyType, valueType), lookupType, _) if keyType.isCoerceableFrom(lookupType) => valueType.validNel
        case (WomCompositeType(typeMap, _), WomStringType, StringLiteral(str)) => typeMap.get(str) match {
          case Some(innerType) => innerType.validNel
          case None => s"Type evaluation failed. No such field '$str' for expression $a".invalidNel
        }
        case (WomObjectType, WomStringType, _) => WomAnyType.validNel
        case (WomAnyType, _, _) => WomAnyType.validNel
        case (otherObject, otherKey, _) => s"Type evaluation failed for $a. Cannot dereference a ${otherObject.stableName} value using a ${otherKey.stableName} key".invalidNel
      }
    }
  }

  /**
    * Try to perform the first lookup in the chain.
    * Then, depending on whether we're done or not, return the result or recurse for the next lookup
    *
    * @param womType The type to perform the lookup on
    * @param lookupChain The chain of lookups to perform
    *
    * @return The ultimate value of the final lookup (or any error along the way!).
    */
  private def doLookup(womType: WomType, lookupChain: NonEmptyList[String]): ErrorOr[WomType] = {
    val key = lookupChain.head
    val tail = NonEmptyList.fromList(lookupChain.tail)

    val thisValue: ErrorOr[WomType] = womType match {
      case WomCompositeType(typeMap, _) => typeMap.get(key).toErrorOr(s"No such field '$key' on type ${womType.stableName}.")
      case WomObjectType => WomAnyType.validNel
      case WomPairType(left, _) if key == "left" => left.validNel
      case WomPairType(_, right) if key == "right" => right.validNel
      case WomMapType(_, right) => right.validNel
      case WomAnyType => WomAnyType.validNel
      case _ => s"No such field '$key' on type ${womType.stableName}".invalidNel
    }

    tail match {
      case None => thisValue
      case Some(lookupTail) => thisValue flatMap { doLookup(_, lookupTail) }
    }
  }
}
