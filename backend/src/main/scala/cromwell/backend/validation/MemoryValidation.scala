package cromwell.backend.validation

import cats.syntax.validated._
import cromwell.backend.MemorySize
import lenthall.validation.ErrorOr._
import wdl4s.parser.MemoryUnit
import wdl4s.types.{WdlIntegerType, WdlStringType}
import wdl4s.values.{WdlInteger, WdlString}

/**
  * Validates the "memory" runtime attribute as an Integer or String with format '8 GB', returning the value as a
  * `MemorySize`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * There is no default, however `optional` can be used to validate the attribute and return the validated value as an
  * `Option`, wrapped in an `Some`, if present, or `None` if not found.
  *
  * `withDefaultMemory` can be used to create a memory validation that defaults to a particular memory size.
  */
object MemoryValidation {
  val key = RuntimeAttributesKeys.MemoryKey

  lazy val instance = new MemoryValidation

  lazy val optional = instance.optional

  def withDefaultMemory(memorySize: MemorySize) = instance.withDefault(WdlInteger(memorySize.bytes.toInt))

  private[validation] val missingMessage =
    s"Expecting $key runtime attribute to be an Integer or String with format '8 GB'"
  private[validation] val wrongAmountFormat = s"Expecting $key runtime attribute value greater than 0 but got %s"
  private[validation] val missingFormat = s"$missingMessage. Exception: %s"

  private[validation] def validateMemoryString(wdlString: WdlString): ErrorOr[MemorySize] =
    validateMemoryString(wdlString.value)

  private[validation] def validateMemoryString(value: String): ErrorOr[MemorySize] = {
    MemorySize.parse(value) match {
      case scala.util.Success(memorySize: MemorySize) if memorySize.amount > 0 =>
        memorySize.to(MemoryUnit.GB).validNel
      case scala.util.Success(memorySize: MemorySize) =>
        wrongAmountFormat.format(memorySize.amount).invalidNel
      case scala.util.Failure(throwable) =>
        missingFormat.format(throwable.getMessage).invalidNel
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

class MemoryValidation extends RuntimeAttributesValidation[MemorySize] {

  import MemoryValidation._

  override def key = RuntimeAttributesKeys.MemoryKey

  override def coercion = Seq(WdlIntegerType, WdlStringType)

  override protected def validateValue = {
    case WdlInteger(value) => MemoryValidation.validateMemoryInteger(value)
    case WdlString(value) => MemoryValidation.validateMemoryString(value)
  }

  override def failureMessage = missingMessage
}
