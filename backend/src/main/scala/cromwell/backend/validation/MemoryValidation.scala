package cromwell.backend.validation

import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr._
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.types.{WomIntegerType, WomLongType, WomStringType}
import wom.values.{WomInteger, WomLong, WomString, WomValue}

import scala.util.{Failure, Success}

/**
  * Validates the "memory" runtime attribute as an Integer or String with format '8 GB', returning the value as a
  * `MemorySize`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * `configDefaultWdlValue` returns the value of the attribute as specified by the
  * reference.conf file, coerced into a WomValue.
  *
  * `optional` can be used to return the validated value as an `Option`,
  * wrapped in a `Some`, if present, or `None` if not found.
  *
  * `withDefaultMemory` can be used to create a memory validation that defaults to a particular memory size.
  */
object MemoryValidation {
  def instance(attributeName: String = RuntimeAttributesKeys.MemoryKey): RuntimeAttributesValidation[MemorySize] =
    new MemoryValidation(attributeName)
  def optional(attributeName: String = RuntimeAttributesKeys.MemoryKey): OptionalRuntimeAttributesValidation[MemorySize] =
    instance(attributeName).optional
  def configDefaultString(attributeName: String = RuntimeAttributesKeys.MemoryKey, config: Option[Config]): Option[String] =
    instance(attributeName).configDefaultValue(config)
  def withDefaultMemory(attributeName: String = RuntimeAttributesKeys.MemoryKey, memorySize: String): RuntimeAttributesValidation[MemorySize] = {
    MemorySize.parse(memorySize) match {
      case Success(memory) => instance(attributeName).withDefault(WomInteger(memory.bytes.toInt))
      case Failure(_) => instance(attributeName).withDefault(BadDefaultAttribute(WomString(memorySize.toString)))
    }
  }

  private[validation] val wrongAmountFormat =
    "Expecting %s runtime attribute value greater than 0 but got %s"
  private[validation] val wrongTypeFormat =
    "Expecting %s runtime attribute to be an Integer or String with format '8 GB'." +
      " Exception: %s"

  private[validation] def validateMemoryString(attributeName: String, wdlString: WomString): ErrorOr[MemorySize] =
    validateMemoryString(attributeName, wdlString.value)

  private[validation] def validateMemoryString(attributeName: String, value: String): ErrorOr[MemorySize] = {
    MemorySize.parse(value) match {
      case scala.util.Success(memorySize: MemorySize) if memorySize.amount > 0 =>
        memorySize.to(MemoryUnit.GB).validNel
      case scala.util.Success(memorySize: MemorySize) =>
        wrongAmountFormat.format(attributeName, memorySize.amount).invalidNel
      case scala.util.Failure(throwable) =>
        wrongTypeFormat.format(attributeName, throwable.getMessage).invalidNel
    }
  }

  private[validation] def validateMemoryInteger(attributeName: String, wdlInteger: WomInteger): ErrorOr[MemorySize] =
    validateMemoryInteger(attributeName, wdlInteger.value)

  private[validation] def validateMemoryInteger(attributeName: String, value: Int): ErrorOr[MemorySize] = {
    if (value <= 0)
      wrongAmountFormat.format(attributeName, value).invalidNel
    else
      MemorySize(value.toDouble, MemoryUnit.Bytes).to(MemoryUnit.GB).validNel
  }

  def validateMemoryLong(attributeName: String, value: Long): ErrorOr[MemorySize] = {
    if (value <= 0)
      wrongAmountFormat.format(attributeName, value).invalidNel
    else
      MemorySize(value.toDouble, MemoryUnit.Bytes).to(MemoryUnit.GB).validNel
  }
}

class MemoryValidation(attributeName: String = RuntimeAttributesKeys.MemoryKey) extends RuntimeAttributesValidation[MemorySize] {

  import MemoryValidation._

  override def key = attributeName

  override def coercion = Seq(WomIntegerType, WomLongType, WomStringType)

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[MemorySize]] = {
    case WomLong(value) => MemoryValidation.validateMemoryLong(key, value)
    case WomInteger(value) => MemoryValidation.validateMemoryInteger(key, value)
    case WomString(value) => MemoryValidation.validateMemoryString(key, value)
  }

  override def missingValueMessage: String = wrongTypeFormat.format(key, "Not supported WDL type value")
}
