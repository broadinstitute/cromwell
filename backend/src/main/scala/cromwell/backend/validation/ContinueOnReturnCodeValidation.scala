package cromwell.backend.validation

import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import com.typesafe.config.Config
import cromwell.backend.validation.RuntimeAttributesValidation._
import lenthall.validation.ErrorOr._
import wdl4s.types.{WdlArrayType, WdlIntegerType, WdlStringType, WdlType}
import wdl4s.values.{WdlArray, WdlBoolean, WdlInteger, WdlString, WdlValue}

import scala.util.Try

/**
  * Validates the "continueOnReturnCode" runtime attribute a Boolean, a String 'true' or 'false', or an Array[Int],
  * returning the value as a instance of a class extending trait `ContinueOnReturnCode`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * `configDefaultWdlValue` returns the value of the attribute as specified by the
  * reference.conf file, coerced into a WdlValue.
  *
  * `default` a validation with the default value specified by the reference.conf file.
  *
  * `optional` can be used to return the validated value as an `Option`,
  * wrapped in a `Some`, if present, or `None` if not found.
  */
object ContinueOnReturnCodeValidation {
  lazy val instance: RuntimeAttributesValidation[ContinueOnReturnCode] = new ContinueOnReturnCodeValidation
  def default(runtimeConfig: Config): RuntimeAttributesValidation[ContinueOnReturnCode] = instance.withDefault(
    configDefaultWdlValue(runtimeConfig) getOrElse WdlInteger(0))
  def configDefaultWdlValue(runtimeConfig: Config): Option[WdlValue] = instance.configDefaultWdlValue(runtimeConfig)
}

class ContinueOnReturnCodeValidation extends RuntimeAttributesValidation[ContinueOnReturnCode] {

  override def key: String = RuntimeAttributesKeys.ContinueOnReturnCodeKey

  override def coercion: Set[WdlType] = ContinueOnReturnCode.validWdlTypes

  override def validateValue: PartialFunction[WdlValue, ErrorOr[ContinueOnReturnCode]] = {
    case WdlBoolean(value) => ContinueOnReturnCodeFlag(value).validNel
    case WdlString(value) if Try(value.toBoolean).isSuccess => ContinueOnReturnCodeFlag(value.toBoolean).validNel
    case WdlString(value) if Try(value.toInt).isSuccess => ContinueOnReturnCodeSet(Set(value.toInt)).validNel
    case WdlInteger(value) => ContinueOnReturnCodeSet(Set(value)).validNel
    case value@WdlArray(_, seq) =>
      val errorOrInts: ErrorOr[List[Int]] = (seq.toList map validateInt).sequence[ErrorOr, Int]
      errorOrInts match {
        case Valid(ints) => ContinueOnReturnCodeSet(ints.toSet).validNel
        case Invalid(_) => invalidValueFailure(value)
      }
  }

  override def validateExpression: PartialFunction[WdlValue, Boolean] = {
    case WdlBoolean(_) => true
    case WdlString(value) if Try(value.toInt).isSuccess => true
    case WdlString(value) if Try(value.toBoolean).isSuccess => true
    case WdlInteger(_) => true
    case WdlArray(WdlArrayType(WdlStringType), elements) =>
      elements forall {
        value => Try(value.valueString.toInt).isSuccess
      }
    case WdlArray(WdlArrayType(WdlIntegerType), _) => true
  }

  override protected def missingValueMessage: String = s"Expecting $key" +
    " runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]"

  override protected def usedInCallCaching: Boolean = true
}
