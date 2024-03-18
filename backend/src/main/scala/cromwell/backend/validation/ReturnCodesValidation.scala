package cromwell.backend.validation

import cats.data.Validated.{Invalid, Valid}
import cats.implicits.{catsSyntaxValidatedId, toTraverseOps}
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.validation.RuntimeAttributesValidation.validateInt
import wom.RuntimeAttributesKeys
import wom.RuntimeAttributesKeys._
import wom.types.{WomArrayType, WomIntegerType, WomStringType, WomType}
import wom.values.{WomArray, WomInteger, WomString, WomValue}

import scala.util.Try

/**
 * Validates the "returnCodes" runtime attribute as an Integer, returning the value as an `Int`, or as
 *  an array of Integers, returning the value as Array[Int], or "*", returning the value as a String.
 *
 * `default` a hardcoded default WomValue for returnCodes.
 *
 * `configDefaultWdlValue` returns the value of the attribute as specified by the
 * reference.conf file, coerced into a WomValue.
 */
object ReturnCodesValidation {
  lazy val instance: RuntimeAttributesValidation[ReturnCodes] = new ReturnCodesValidation(ReturnCodesKey)
  lazy val optional: OptionalRuntimeAttributesValidation[ReturnCodes] = instance.optional
  def default(runtimeConfig: Option[Config]): RuntimeAttributesValidation[ReturnCodes] =
    instance.withDefault(configDefaultWdlValue(runtimeConfig) getOrElse WomInteger(0))

  def configDefaultWdlValue(runtimeConfig: Option[Config]): Option[WomValue] =
    instance.configDefaultWomValue(runtimeConfig)
  def configDefaultWomValue(config: Option[Config]): Option[WomValue] = instance.configDefaultWomValue(config)
}

class ReturnCodesValidation(attributeName: String) extends RuntimeAttributesValidation[ReturnCodes] {

  override def key: String = RuntimeAttributesKeys.ReturnCodesKey

  override def coercion: Set[WomType] = ReturnCodes.validWdlTypes

  override def validateValue: PartialFunction[WomValue, ErrorOr[ReturnCodes]] = {
    case WomString(value) if value.equals("*") =>
      ReturnCodesString(value).validNel
    case WomString(value) if Try(value.toInt).isSuccess =>
      ReturnCodesSet(Set(value.toInt)).validNel
    case WomInteger(value) =>
      ReturnCodesSet(Set(value)).validNel
    case value @ WomArray(_, seq) =>
      val errorOrInts: ErrorOr[List[Int]] = (seq.toList map validateInt).sequence[ErrorOr, Int]
      errorOrInts match {
        case Valid(ints) => ReturnCodesSet(ints.toSet).validNel
        case Invalid(_) => invalidValueFailure(value)
      }
  }

  override def validateExpression: PartialFunction[WomValue, Boolean] = {
    case WomString(value) if Try(value.toInt).isSuccess => true
    case WomInteger(_) => true
    case WomString(value) if value.equals("*") => true
    case WomArray(WomArrayType(WomStringType), elements) =>
      elements forall { value =>
        Try(value.valueString.toInt).isSuccess
      }
    case WomArray(WomArrayType(WomIntegerType), _) => true
  }

  override protected def missingValueMessage: String = s"Expecting $key" +
    " runtime attribute to be either a String '*' or an Array[Int]"

  override def usedInCallCaching: Boolean = true
}
