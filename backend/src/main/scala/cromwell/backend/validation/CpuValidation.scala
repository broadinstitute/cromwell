package cromwell.backend.validation

import cats.syntax.validated._
import com.typesafe.config.Config
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wdl.types.WdlIntegerType
import wdl4s.wdl.values.{WdlInteger, WdlValue}

/**
  * Validates the "cpu" runtime attribute an Integer greater than 0, returning the value as an `Int`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * `default` a hardcoded default WdlValue for Cpu.
  *
  * `configDefaultWdlValue` returns the value of the attribute as specified by the
  * reference.conf file, coerced into a WdlValue.
  */
object CpuValidation {
  lazy val instance: RuntimeAttributesValidation[Int] = new CpuValidation
  lazy val optional: OptionalRuntimeAttributesValidation[Int] = instance.optional
  lazy val default: WdlValue = WdlInteger(1)
  def configDefaultWdlValue(config: Option[Config]): Option[WdlValue] = instance.configDefaultWdlValue(config)
}

class CpuValidation extends IntRuntimeAttributesValidation(RuntimeAttributesKeys.CpuKey) {
  override protected def validateValue: PartialFunction[WdlValue, ErrorOr[Int]] = {
    case wdlValue if WdlIntegerType.coerceRawValue(wdlValue).isSuccess =>
      WdlIntegerType.coerceRawValue(wdlValue).get match {
        case WdlInteger(value) =>
          if (value.toInt <= 0)
            s"Expecting $key runtime attribute value greater than 0".invalidNel
          else
            value.toInt.validNel
      }
  }

  override protected def missingValueMessage: String = s"Expecting $key runtime attribute to be an Integer"
}
