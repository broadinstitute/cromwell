package wdl.transforms.base.linking.expression.values

import cats.data.{NonEmptyList, Validated}
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{EvaluatedValue, ForCommandInstantiationOptions, ValueEvaluator}
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.types._
import wom.values._
import wdl.transforms.base.wdlom2wdl.WdlWriter.ops._
import wdl.transforms.base.wdlom2wdl.WdlWriterImpl._


object LookupEvaluators {

  implicit val identifierLookupEvaluator: ValueEvaluator[IdentifierLookup] = new ValueEvaluator[IdentifierLookup] {
    override def evaluateValue(a: IdentifierLookup,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
      inputs.get(a.identifier) match {
        case Some(value) =>
          val mapped = forCommandInstantiationOptions.fold(value)(_.valueMapper(value))
          EvaluatedValue(mapped, Seq.empty).validNel
        case None =>
          s"ValueEvaluator[IdentifierLookup]: No suitable input for '${a.identifier}' amongst {${inputs.keys.mkString(", ")}}".invalidNel
      }
    }
  }

  implicit val expressionMemberAccessEvaluator: ValueEvaluator[ExpressionMemberAccess] = new ValueEvaluator[ExpressionMemberAccess] {
    override def evaluateValue(a: ExpressionMemberAccess,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
      a.expression.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions) flatMap { evaluated =>
        doLookup(evaluated.value, a.memberAccessTail) map { EvaluatedValue(_, evaluated.sideEffectFiles) }
      }
    }
  }

  implicit val identifierMemberAccessEvaluator: ValueEvaluator[IdentifierMemberAccess] = new ValueEvaluator[IdentifierMemberAccess] {
    override def evaluateValue(a: IdentifierMemberAccess,
                               inputs: Map[String, WomValue],
                               ioFunctionSet: IoFunctionSet,
                               forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {

      // Do the first lookup and decide whether any more lookups are needed:
      val generatedValueAndLookups: ErrorOr[(WomValue, Seq[String])] = {
        val callOutputKey = s"${a.first}.${a.second}"

        if (inputs.keySet.contains(a.first)) { (inputs(a.first), List(a.second) ++ a.memberAccessTail).validNel }
        else if (inputs.keySet.contains(callOutputKey)) { (inputs(callOutputKey), a.memberAccessTail).validNel }
        else
          s"No value found for member access lookup. Report this bug: Insufficient input values supplied by engine. Needed '${a.first}' or '$callOutputKey' but only received: '${inputs.keys.mkString(", ")}'".invalidNel
      }

      generatedValueAndLookups flatMap { case (foundValue, lookups) => NonEmptyList.fromList(lookups.toList) match {
        case Some(lookupNel) => doLookup(foundValue, lookupNel) map { EvaluatedValue(_, Seq.empty) }
        case None => EvaluatedValue(foundValue, Seq.empty).validNel
      }}
    }
  }

  implicit val indexAccessValueEvaluator: ValueEvaluator[IndexAccess] = new ValueEvaluator[IndexAccess] {
    override def evaluateValue(a: IndexAccess, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, forCommandInstantiationOptions: Option[ForCommandInstantiationOptions])
                              (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[EvaluatedValue[_ <: WomValue]] = {
      (a.expressionElement.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions), a.index.evaluateValue(inputs, ioFunctionSet, forCommandInstantiationOptions)) flatMapN { (lhs, rhs) =>
        val value: ErrorOr[WomValue] = (lhs.value, rhs.value) match {
          case (array: WomArray, WomInteger(index)) =>
            if (array.value.length > index)
              array.value(index).validNel
            else
              s"Bad array access ${a.toWdlV1}: Array size ${array.value.length} does not have an index value '$index'".invalidNel
          case (WomObject(values, _), WomString(index)) =>
            if(values.contains(index))
              values(index).validNel
            else
              s"Bad Object access ${a.toWdlV1}: Object with keys [${values.keySet.mkString(", ")}] does not have an index value [$index]".invalidNel
          case (WomMap(mapType, values), index) =>
            if(values.contains(index))
              values(index).validNel
            else
              s"Bad Map access ${a.toWdlV1}: This ${mapType.stableName} does not have a ${index.womType.stableName} index value [${index.toWomString}]".invalidNel
          case (otherCollection, otherKey) =>
            s"Bad index access ${a.toWdlV1}: Cannot use '${otherKey.womType.stableName}' to index '${otherCollection.womType.stableName}'".invalidNel
        }

        value map { EvaluatedValue(_, lhs.sideEffectFiles ++ rhs.sideEffectFiles) }
      }
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
      case WomObject(_, WomCompositeType(typeMap, _)) if typeMap.contains(key) => typeMap(key) match {
        case WomOptionalType(innerType) => WomOptionalValue(innerType, None).validNel
        case other => s"Composite value was unexpectedly missing a field: '$key' (expected type ${other.stableName}). Report this bug! Static validation failed.".invalidNel
      }
      case WomObject(_, _: WomCompositeType) => s"No such field '$key' on type ${womValue.womType.stableName}. Report this bug! Static validation failed.".invalidNel
      case WomObject(_, _) => s"'Object'-type value did not contain the field '$key' at runtime".invalidNel
      case p: WomPair if key == "left" => p.left.validNel
      case p: WomPair if key == "right" => p.right.validNel
      case WomMap(_, value) => Validated.fromOption(
        o = value.collectFirst {
          case (k, v) if k.valueString == key => v
        },
        ifNone = NonEmptyList(s"Requested key '$key' not found in the Map. Available keys were: ${value.keySet.mkString("[ ", ",", "]")}", Nil))
      case _ => s"No such field '$key' on type ${womValue.womType.stableName}. Report this bug! Static validation failed.".invalidNel
    }

    tail match {
      case None => thisValue
      case Some(lookupTail) => thisValue flatMap { doLookup(_, lookupTail) }
    }

  }

}
