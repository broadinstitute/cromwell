package cromwell.backend.validation

import cats.data.Validated.{Invalid, Valid}
import cats.implicits.{catsSyntaxValidatedId, toTraverseOps}
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.validation.RuntimeAttributesValidation.validateInt
import wom.types.{WomArrayType, WomIntegerType, WomStringType, WomType}
import wom.values.{WomArray, WomInteger, WomString, WomValue}

import scala.util.Try

/**
 * Validates the `continueOnReturnCode` and `returnCodes` runtime attributes, since they share the functionality of
 * being a Set of Integers. 
 */
trait ReturnCodeValidation extends RuntimeAttributesValidation[ReturnCode] {
  override def key: String
  override def coercion: Iterable[WomType]
  override def validateValue: PartialFunction[WomValue, ErrorOr[ReturnCode]] = {
    case WomString(value) if Try(value.toInt).isSuccess => ReturnCodeSet(Set(value.toInt)).validNel
    case WomInteger(value) => ReturnCodeSet(Set(value)).validNel
    case value @ WomArray(_, seq) =>
      val errorOrInts: ErrorOr[List[Int]] = (seq.toList map validateInt).sequence[ErrorOr, Int]
      errorOrInts match {
        case Valid(ints) => ReturnCodeSet(ints.toSet).validNel
        case Invalid(_) => invalidValueFailure(value)
      }
    case value => invalidValueFailure(value)
  }

  override def validateExpression: PartialFunction[WomValue, Boolean] = {
    case WomString(value) if Try(value.toInt).isSuccess => true
    case WomInteger(_) => true
    case WomArray(WomArrayType(WomStringType), elements) =>
      elements forall { value =>
        Try(value.valueString.toInt).isSuccess
      }
    case WomArray(WomArrayType(WomIntegerType), _) => true
    case _ => false
  }

  override def usedInCallCaching: Boolean = true
}
