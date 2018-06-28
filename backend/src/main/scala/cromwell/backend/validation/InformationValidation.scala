package cromwell.backend.validation

import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr._
import squants.QuantityParseException
import squants.information.{Information, InformationUnit}
import wom.types.{WomIntegerType, WomLongType, WomStringType}
import wom.values.{WomInteger, WomLong, WomString, WomValue}

import scala.util.{Failure, Success}

/**
  * Validates a runtime attribute expressed as information (bytes)
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * `configDefaultWdlValue` returns the value of the attribute as specified by the
  * reference.conf file, coerced into a WomValue.
  *
  * `optional` can be used to return the validated value as an `Option`,
  * wrapped in a `Some`, if present, or `None` if not found.
  *
  * `withDefault` can be used to create a validation that defaults to a particular size.
  */
object InformationValidation {
  def instance(attributeName: String, defaultUnit: InformationUnit): RuntimeAttributesValidation[Information] =
    new InformationValidation(attributeName, defaultUnit)
  def optional(attributeName: String, defaultUnit: InformationUnit): OptionalRuntimeAttributesValidation[Information] =
    instance(attributeName, defaultUnit).optional
  def configDefaultString(attributeName: String, config: Option[Config], defaultUnit: InformationUnit): Option[String] =
    instance(attributeName, defaultUnit).configDefaultValue(config)
  def withDefault(attributeName: String, default: String, defaultUnit: InformationUnit): RuntimeAttributesValidation[Information] = {
    Information(default) match {
      case Success(information) => instance(attributeName, defaultUnit).withDefault(WomInteger(information.toBytes.toInt))
      case Failure(_) => instance(attributeName, defaultUnit).withDefault(BadDefaultAttribute(WomString(default.toString)))
    }
  }

  private[validation] val wrongAmountFormat =
    "Expecting %s runtime attribute value greater than 0 but got %s"
  private[validation] val wrongTypeFormat =
    "Expecting %s runtime attribute to be an Integer or String with format '8 GB'." +
      " Exception: %s"

  private[validation] def validateString(attributeName: String, wdlString: WomString): ErrorOr[Information] =
    validateString(attributeName, wdlString.value)

  private[validation] def validateString(attributeName: String, value: String): ErrorOr[Information] = {
    Information(value) match {
      case scala.util.Success(information: Information) if information.value > 0D =>
        information.validNel
      case scala.util.Success(information: Information) =>
        wrongAmountFormat.format(attributeName, information.value).invalidNel
      case scala.util.Failure(_: QuantityParseException) =>
        wrongTypeFormat.format(attributeName, s"$value should be of the form 'X Unit' where X is a number, e.g. 8 GB").invalidNel
      case scala.util.Failure(throwable) =>
        wrongTypeFormat.format(attributeName, throwable.getMessage).invalidNel
    }
  }

  private[validation] def validateInteger(attributeName: String, value: Int, defaultUnit: InformationUnit): ErrorOr[Information] = {
    if (value <= 0)
      wrongAmountFormat.format(attributeName, value).invalidNel
    else
      defaultUnit(value.toDouble).validNel
  }

  def validateLong(attributeName: String, value: Long, defaultUnit: InformationUnit): ErrorOr[Information] = {
    if (value <= 0)
      wrongAmountFormat.format(attributeName, value).invalidNel
    else
      defaultUnit(value.toDouble).validNel
  }
}

class InformationValidation(attributeName: String, defaultUnit: InformationUnit) extends RuntimeAttributesValidation[Information] {

  import InformationValidation._

  override def key = attributeName

  override def coercion = Seq(WomIntegerType, WomLongType, WomStringType)

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Information]] = {
    case WomLong(value) => InformationValidation.validateLong(key, value, defaultUnit)
    case WomInteger(value) => InformationValidation.validateInteger(key, value, defaultUnit)
    case WomString(value) => InformationValidation.validateString(key, value)
  }

  override def missingValueMessage: String = wrongTypeFormat.format(key, "Not supported WDL type value")
}
