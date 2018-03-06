package wdl.draft3.transforms.linking.expression.values

import cats.data.NonEmptyList
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.ValueEvaluator
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.model.draft3.graph._
import wom.expression.IoFunctionSet
import wom.types.{WomCompositeType, WomOptionalType}
import wom.values.{WomObject, WomOptionalValue, WomPair, WomValue}


object LookupEvaluators {

  implicit val identifierLookupEvaluator: ValueEvaluator[IdentifierLookup] = new ValueEvaluator[IdentifierLookup] {
    override def evaluateValue(a: IdentifierLookup,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      inputs.get(a.identifier) match {
        case Some(value) => value.validNel
        case None => s"No suitable input for identifier lookup '${a.identifier}' amongst {${inputs.keys.mkString(", ")}}".invalidNel
      }
    }
  }

  implicit val expressionMemberAccessEvaluator: ValueEvaluator[ExpressionMemberAccess] = new ValueEvaluator[ExpressionMemberAccess] {
    override def evaluateValue(a: ExpressionMemberAccess,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      a.expression.evaluateValue(inputs, ioFunctionSet) flatMap { doLookup(_, a.memberAccessTail) }
    }
  }

  implicit val identifierMemberAccessEvaluator: ValueEvaluator[IdentifierMemberAccess] = new ValueEvaluator[IdentifierMemberAccess] {
    override def evaluateValue(a: IdentifierMemberAccess,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {

      // Do the first lookup and decide whether any more lookups are needed:
      val generatedValueAndLookups: ErrorOr[(WomValue, Seq[String])] = {
        val callOutputKey = s"${a.first}.${a.second}"

        if (inputs.keySet.contains(a.first)) { (inputs(a.first), List(a.second) ++ a.memberAccessTail).validNel }
        else if (inputs.keySet.contains(callOutputKey)) { (inputs(callOutputKey), a.memberAccessTail).validNel }
        else
          s"No value found for member access lookup. Report this bug: Insufficient input values supplied by engine. Needed '${a.first}' or '$callOutputKey' but only received: '${inputs.keys.mkString(", ")}'".invalidNel
      }

      generatedValueAndLookups flatMap { case (foundValue, lookups) => NonEmptyList.fromList(lookups.toList) match {
        case Some(lookupNel) => doLookup(foundValue, lookupNel)
        case None => foundValue.validNel
      }}
    }
  }

  /**
    * Try to perform the first lookup in the chain.
    * Then, depending on whether we're done or not, return the result or recurse for the next lookup
    *
    * @param womValue The value to perform the lookup on
    * @param lookupChain The chain of lookups to perform
    *
    * @return The ultimate value of the final lookup (or any error along the way!).
    */
  private def doLookup(womValue: WomValue, lookupChain: NonEmptyList[String]): ErrorOr[WomValue] = {
    val key = lookupChain.head
    val tail = NonEmptyList.fromList(lookupChain.tail)

    val thisValue: ErrorOr[WomValue] = womValue match {
      case WomObject(values, _) if values.contains(key) => values(key).validNel
      case WomObject(_, WomCompositeType(typeMap)) if typeMap.contains(key) => typeMap(key) match {
        case WomOptionalType(innerType) => WomOptionalValue(innerType, None).validNel
        case other => s"Composite value was unexpectedly missing a field: '$key' (expected type ${other.toDisplayString}). Report this bug! Static validation failed.".invalidNel
      }
      case WomObject(_, _: WomCompositeType) => s"No such field '$key' on type ${womValue.womType.toDisplayString}. Report this bug! Static validation failed.".invalidNel
      case WomObject(_, _) => s"'Object'-type value did not contain the field '$key' at runtime".invalidNel
      case p: WomPair if key == "left" => p.left.validNel
      case p: WomPair if key == "right" => p.right.validNel
      case _ => s"No such field '$key' on type ${womValue.womType.toDisplayString}. Report this bug! Static validation failed.".invalidNel
    }

    tail match {
      case None => thisValue
      case Some(lookupTail) => thisValue flatMap { doLookup(_, lookupTail) }
    }

  }

}
