package cromwell.backend.validation

import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
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
  lazy val instance: RuntimeAttributesValidation[Int] = new CpuValidation
  lazy val optional: OptionalRuntimeAttributesValidation[Int] = instance.optional
  lazy val default: WomValue = WomInteger(1)
  def configDefaultWdlValue(config: Option[Config]): Option[WomValue] = instance.configDefaultWdlValue(config)
}

class CpuValidation extends IntRuntimeAttributesValidation(RuntimeAttributesKeys.CpuKey) {
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
