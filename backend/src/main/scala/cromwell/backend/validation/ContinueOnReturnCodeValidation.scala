package cromwell.backend.validation

import cats.implicits._
import com.typesafe.config.Config
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
  lazy val instance: RuntimeAttributesValidation[ReturnCode] = new ContinueOnReturnCodeValidation
  def default(runtimeConfig: Option[Config]): RuntimeAttributesValidation[ReturnCode] =
    instance.withDefault(configDefaultWdlValue(runtimeConfig) getOrElse WomInteger(0))
  def configDefaultWdlValue(runtimeConfig: Option[Config]): Option[WomValue] =
    instance.configDefaultWomValue(runtimeConfig)
}

class ContinueOnReturnCodeValidation extends ReturnCodeValidation {

  override def key: String = RuntimeAttributesKeys.ContinueOnReturnCodeKey

  override def coercion: Set[WomType] = ContinueOnReturnCode.validWdlTypes

  override def validateValue: PartialFunction[WomValue, ErrorOr[ReturnCode]] = {
    case WomBoolean(value) => ContinueOnReturnCodeFlag(value).validNel
    case WomString(value) if Try(value.toBoolean).isSuccess => ContinueOnReturnCodeFlag(value.toBoolean).validNel
    case value => super.validateValue(value)
  }

  override def validateExpression: PartialFunction[WomValue, Boolean] = {
    case WomBoolean(_) => true
    case WomString(value) if Try(value.toBoolean).isSuccess => true
    case value => super.validateExpression(value)
  }

  override protected def missingValueMessage: String = s"Expecting $key" +
    " runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]"

  override def usedInCallCaching: Boolean = true
}
