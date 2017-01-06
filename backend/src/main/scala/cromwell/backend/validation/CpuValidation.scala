package cromwell.backend.validation

import cats.syntax.validated._
import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.types.WdlIntegerType
import wdl4s.values.{WdlInteger, WdlValue}

/**
  * Validates the "cpu" runtime attribute an Integer greater than 0, returning the value as an `Int`.
  *
  * `instance` returns an validation that errors when no attribute is specified.
  *
  * The default returns `1` when no attribute is specified.
  *
  * `optional` can be used return the validated value as an `Option`, wrapped in a `Some`, if present, or `None` if not
  * found.
  */
object CpuValidation extends {
  lazy val instance: RuntimeAttributesValidation[Int] = new CpuValidation
  lazy val default: RuntimeAttributesValidation[Int] = instance.withDefault(WdlInteger(1))
  lazy val optional: OptionalRuntimeAttributesValidation[Int] = default.optional
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
