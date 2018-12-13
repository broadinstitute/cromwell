package cromwell.backend.validation

import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr._
import wdl4s.parser.MemoryUnit
import wom.RuntimeAttributesKeys
import wom.format.MemorySize
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
  def instance(attributeName: String = RuntimeAttributesKeys.MemoryKey, defaultUnit: MemoryUnit, allowZero: Boolean = false): RuntimeAttributesValidation[MemorySize] =
    new InformationValidation(attributeName, defaultUnit, allowZero)
  def optional(attributeName: String = RuntimeAttributesKeys.MemoryKey, defaultUnit: MemoryUnit, allowZero: Boolean = false): OptionalRuntimeAttributesValidation[MemorySize] =
    instance(attributeName, defaultUnit, allowZero).optional
  def configDefaultString(attributeName: String = RuntimeAttributesKeys.MemoryKey, config: Option[Config], defaultUnit: MemoryUnit, allowZero: Boolean = false): Option[String] =
    instance(attributeName, defaultUnit, allowZero).configDefaultValue(config)
  def withDefaultMemory(attributeName: String = RuntimeAttributesKeys.MemoryKey, memorySize: String, defaultUnit: MemoryUnit, allowZero: Boolean = false): RuntimeAttributesValidation[MemorySize] = {
    MemorySize.parse(memorySize) match {
      case Success(memory) => instance(attributeName, defaultUnit, allowZero).withDefault(WomLong(memory.bytes.toLong))
      case Failure(_) => instance(attributeName, defaultUnit, allowZero).withDefault(BadDefaultAttribute(WomString(memorySize.toString)))
    }
  }

  private[validation] val wrongAmountFormat =
    "Expecting %s runtime attribute value greater than 0 but got %s"
  private[validation] val wrongTypeFormat =
    "Expecting %s runtime attribute to be an Integer or String with format '8 GB'." +
      " Exception: %s"

  private[validation] def validateString(attributeName: String, wdlString: WomString, allowZero: Boolean): ErrorOr[MemorySize] =
    validateString(attributeName, wdlString.value, allowZero)

  private[validation] def validateString(attributeName: String, value: String, allowZero: Boolean): ErrorOr[MemorySize] = {
    MemorySize.parse(value) match {
      case scala.util.Success(memorySize: MemorySize) if memorySize.amount > 0 || (memorySize.amount == 0 && allowZero) =>
        memorySize.to(MemoryUnit.GB).validNel
      case scala.util.Success(memorySize: MemorySize) =>
        wrongAmountFormat.format(attributeName, memorySize.amount).invalidNel
      case scala.util.Failure(throwable) =>
        wrongTypeFormat.format(attributeName, throwable.getMessage).invalidNel
    }
  }

  private[validation] def validateInteger(attributeName: String, wdlInteger: WomInteger, defaultUnit: MemoryUnit, allowZero: Boolean): ErrorOr[MemorySize] =
    validateInteger(attributeName, wdlInteger.value, defaultUnit, allowZero)

  private[validation] def validateInteger(attributeName: String, value: Int, defaultUnit: MemoryUnit, allowZero: Boolean): ErrorOr[MemorySize] = {
    if (value < 0 || (value == 0 && !allowZero))
      wrongAmountFormat.format(attributeName, value).invalidNel
    else
      MemorySize(value.toDouble, defaultUnit).to(MemoryUnit.GB).validNel
  }

  def validateLong(attributeName: String, value: Long, defaultUnit: MemoryUnit, allowZero: Boolean): ErrorOr[MemorySize] = {
    if (value < 0 || (value == 0 && !allowZero))
      wrongAmountFormat.format(attributeName, value).invalidNel
    else
      MemorySize(value.toDouble, defaultUnit).to(MemoryUnit.GB).validNel
  }
}

class InformationValidation(attributeName: String, defaultUnit: MemoryUnit, allowZero: Boolean = false) extends RuntimeAttributesValidation[MemorySize] {

  import InformationValidation._

  override def key = attributeName

  override def coercion = Seq(WomIntegerType, WomLongType, WomStringType)

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[MemorySize]] = {
    case WomLong(value) => InformationValidation.validateLong(key, value, defaultUnit, allowZero)
    case WomInteger(value) => InformationValidation.validateInteger(key, value, defaultUnit, allowZero)
    case WomString(value) => InformationValidation.validateString(key, value, allowZero)
  }

  override def missingValueMessage: String = wrongTypeFormat.format(key, "Not supported WDL type value")
}
