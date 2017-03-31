package cromwell.backend.validation

import cats.syntax.validated._
import com.typesafe.config.Config
import cromwell.backend.MemorySize
import lenthall.validation.ErrorOr._
import wdl4s.parser.MemoryUnit
import wdl4s.types.{WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlInteger, WdlString, WdlValue}

import scala.util.{Failure, Success}

/**
  * Validates the "memory" runtime attribute as an Integer or String with format '8 GB', returning the value as a
  * `MemorySize`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * `configDefaultWdlValue` returns the value of the attribute as specified by the
  * reference.conf file, coerced into a WdlValue.
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
      case Success(memory) => instance(attributeName).withDefault(WdlInteger(memory.bytes.toInt))
      case Failure(_) => instance(attributeName).withDefault(BadDefaultAttribute(WdlString(memorySize.toString)))
    }
  }

  private[validation] val wrongAmountFormat =
    s"Expecting ${RuntimeAttributesKeys.MemoryKey} runtime attribute value greater than 0 but got %s"
  private[validation] val wrongTypeFormat =
    s"Expecting ${RuntimeAttributesKeys.MemoryKey} runtime attribute to be an Integer or String with format '8 GB'." +
      s" Exception: %s"

  private[validation] def validateMemoryString(wdlString: WdlString): ErrorOr[MemorySize] =
    validateMemoryString(wdlString.value)

  private[validation] def validateMemoryString(value: String): ErrorOr[MemorySize] = {
    MemorySize.parse(value) match {
      case scala.util.Success(memorySize: MemorySize) if memorySize.amount > 0 =>
        memorySize.to(MemoryUnit.GB).validNel
      case scala.util.Success(memorySize: MemorySize) =>
        wrongAmountFormat.format(memorySize.amount).invalidNel
      case scala.util.Failure(throwable) =>
        wrongTypeFormat.format(throwable.getMessage).invalidNel
    }
  }

  private[validation] def validateMemoryInteger(wdlInteger: WdlInteger): ErrorOr[MemorySize] =
    validateMemoryInteger(wdlInteger.value)

  private[validation] def validateMemoryInteger(value: Int): ErrorOr[MemorySize] = {
    if (value <= 0)
      wrongAmountFormat.format(value).invalidNel
    else
      MemorySize(value.toDouble, MemoryUnit.Bytes).to(MemoryUnit.GB).validNel
  }
}

class MemoryValidation(attributeName: String = RuntimeAttributesKeys.MemoryKey) extends RuntimeAttributesValidation[MemorySize] {

  import MemoryValidation._

  override def key = attributeName

  override def coercion = Seq(WdlIntegerType, WdlStringType)

  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[MemorySize]] = {
    case WdlInteger(value) => MemoryValidation.validateMemoryInteger(value)
    case WdlString(value) => MemoryValidation.validateMemoryString(value)
  }

  override def missingValueMessage: String = wrongTypeFormat.format("Not supported WDL type value")
}
