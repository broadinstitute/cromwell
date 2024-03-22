package cromwell.backend.validation

import cats.implicits.catsSyntaxValidatedId
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import wom.RuntimeAttributesKeys
import wom.types.WomType
import wom.values.{WomInteger, WomString, WomValue}

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
  lazy val instance: RuntimeAttributesValidation[ReturnCode] = new ReturnCodesValidation
  lazy val optional: OptionalRuntimeAttributesValidation[ReturnCode] = instance.optional
  def default(runtimeConfig: Option[Config]): RuntimeAttributesValidation[ReturnCode] =
    instance.withDefault(configDefaultWdlValue(runtimeConfig) getOrElse WomInteger(0))

  def configDefaultWdlValue(runtimeConfig: Option[Config]): Option[WomValue] =
    instance.configDefaultWomValue(runtimeConfig)
  def configDefaultWomValue(config: Option[Config]): Option[WomValue] = instance.configDefaultWomValue(config)
}

class ReturnCodesValidation extends ReturnCodeValidation {

  override def key: String = RuntimeAttributesKeys.ReturnCodesKey

  override def coercion: Set[WomType] = ReturnCodes.validWdlTypes

  override def validateValue: PartialFunction[WomValue, ErrorOr[ReturnCode]] = {
    case WomString(value) if value.equals("*") =>
      ReturnCodesString(value).validNel
    case value => super.validateValue(value)
  }

  override def validateExpression: PartialFunction[WomValue, Boolean] = {
    case WomString(value) if value.equals("*") => true
    case value => super.validateExpression(value)
  }

  override protected def missingValueMessage: String = s"Expecting $key" +
    " runtime attribute to be either a String '*' or an Array[Int]"

  override def usedInCallCaching: Boolean = true
}
