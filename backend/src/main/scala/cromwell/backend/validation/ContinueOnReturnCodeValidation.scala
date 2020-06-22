package cromwell.backend.validation

import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.typesafe.config.Config
import cromwell.backend.validation.RuntimeAttributesValidation._
import common.validation.ErrorOr._
import wom.RuntimeAttributesKeys
import wom.types._
import wom.values._

import scala.util.Try

/**
  * Validates the "continueOnReturnCode" runtime attribute a Boolean, a String 'true' or 'false', or an Array[Int],
  * returning the value as a instance of a class extending trait `ContinueOnReturnCode`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * `configDefaultWdlValue` returns the value of the attribute as specified by the
  * reference.conf file, coerced into a WomValue.
  *
  * `default` a validation with the default value specified by the reference.conf file.
  */
object ContinueOnReturnCodeValidation {
  lazy val instance: RuntimeAttributesValidation[ContinueOnReturnCode] = new ContinueOnReturnCodeValidation
  def default(runtimeConfig: Option[Config]): RuntimeAttributesValidation[ContinueOnReturnCode] = instance.withDefault(
    configDefaultWdlValue(runtimeConfig) getOrElse WomInteger(0))
  def configDefaultWdlValue(runtimeConfig: Option[Config]): Option[WomValue] = instance.configDefaultWomValue(runtimeConfig)
}

class ContinueOnReturnCodeValidation extends RuntimeAttributesValidation[ContinueOnReturnCode] {

  override def key: String = RuntimeAttributesKeys.ContinueOnReturnCodeKey

  override def coercion: Set[WomType] = ContinueOnReturnCode.validWdlTypes

  override def validateValue: PartialFunction[WomValue, ErrorOr[ContinueOnReturnCode]] = {
    case WomBoolean(value) => ContinueOnReturnCodeFlag(value).validNel
    case WomString(value) if Try(value.toBoolean).isSuccess => ContinueOnReturnCodeFlag(value.toBoolean).validNel
    case WomString(value) if Try(value.toInt).isSuccess => ContinueOnReturnCodeSet(Set(value.toInt)).validNel
    case WomInteger(value) => ContinueOnReturnCodeSet(Set(value)).validNel
    case value@WomArray(_, seq) =>
      val errorOrInts: ErrorOr[List[Int]] = (seq.toList map validateInt).sequence[ErrorOr, Int]
      errorOrInts match {
        case Valid(ints) => ContinueOnReturnCodeSet(ints.toSet).validNel
        case Invalid(_) => invalidValueFailure(value)
      }
  }

  override def validateExpression: PartialFunction[WomValue, Boolean] = {
    case WomBoolean(_) => true
    case WomString(value) if Try(value.toInt).isSuccess => true
    case WomString(value) if Try(value.toBoolean).isSuccess => true
    case WomInteger(_) => true
    case WomArray(WomArrayType(WomStringType), elements) =>
      elements forall {
        value => Try(value.valueString.toInt).isSuccess
      }
    case WomArray(WomArrayType(WomIntegerType), _) => true
  }

  override protected def missingValueMessage: String = s"Expecting $key" +
    " runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]"

  override def usedInCallCaching: Boolean = true
}
