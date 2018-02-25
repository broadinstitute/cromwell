package cromwell.backend.validation

import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import wom.RuntimeAttributesKeys._
import wom.types.WomIntegerType
import wom.values.{WomInteger, WomValue}

/**
  * Validates the "cpu" runtime attribute an Integer greater than 0, returning the value as an `Int`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * `default` a hardcoded default WomValue for Cpu.
  *
  * `configDefaultWdlValue` returns the value of the attribute as specified by the
  * reference.conf file, coerced into a WomValue.
  */
object CpuValidation {
  lazy val instance: RuntimeAttributesValidation[Int] = new CpuValidation(CpuKey)
  lazy val optional: OptionalRuntimeAttributesValidation[Int] = instance.optional
  lazy val instanceMin: RuntimeAttributesValidation[Int] = new CpuValidation(CpuMinKey)
  lazy val optionalMin: OptionalRuntimeAttributesValidation[Int] = instanceMin.optional
  lazy val instanceMax: RuntimeAttributesValidation[Int] = new CpuValidation(CpuMaxKey)
  lazy val optionalMax: OptionalRuntimeAttributesValidation[Int] = instanceMax.optional

  lazy val defaultMin: WomValue = WomInteger(1)
  def configDefaultWomValue(config: Option[Config]): Option[WomValue] = instance.configDefaultWomValue(config)
}

class CpuValidation(attributeName: String) extends IntRuntimeAttributesValidation(attributeName) {
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Int]] = {
    case womValue if WomIntegerType.coerceRawValue(womValue).isSuccess =>
      WomIntegerType.coerceRawValue(womValue).get match {
        case WomInteger(value) =>
          if (value.toInt <= 0)
            s"Expecting $key runtime attribute value greater than 0".invalidNel
          else
            value.toInt.validNel
      }
  }

  override protected def missingValueMessage: String = s"Expecting $key runtime attribute to be an Integer"
}
