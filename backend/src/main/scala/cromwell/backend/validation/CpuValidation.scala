package cromwell.backend.validation

import cats.syntax.validated._
import com.typesafe.config.Config
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.types.WdlIntegerType
import wdl4s.values.{WdlInteger, WdlValue}

/**
  * Validates the "cpu" runtime attribute an Integer greater than 0, returning the value as an `Int`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * `configDefaultWdlValue` returns the value of the attribute as specified by the
  * reference.conf file, coerced into a WdlValue.
  *
  * `default` a validation with the default value specified by the reference.conf file.
  *
  * `optional` can be used to return the validated value as an `Option`,
  * wrapped in a `Some`, if present, or `None` if not found.
  */
object CpuValidation extends {
  lazy val instance: RuntimeAttributesValidation[Int] = new CpuValidation
  lazy val default: WdlValue = WdlInteger(1)
  def configDefaultWdlValue(config: Config): Option[WdlValue] = instance.configDefaultWdlValue(config)
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
