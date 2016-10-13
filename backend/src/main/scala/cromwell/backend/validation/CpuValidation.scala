package cromwell.backend.validation

import cats.syntax.validated._
import wdl4s.types.WdlIntegerType
import wdl4s.values.WdlInteger

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
  val key = RuntimeAttributesKeys.CpuKey

  lazy val instance = new CpuValidation

  lazy val default = instance.withDefault(WdlInteger(1))

  lazy val optional = default.optional

  private[validation] val missingMessage = s"Expecting $key runtime attribute to be an Integer"
  private[validation] val wrongAmountMsg = s"Expecting $key runtime attribute value greater than 0"
}

class CpuValidation extends RuntimeAttributesValidation[Int] {

  import CpuValidation._

  override def key = RuntimeAttributesKeys.CpuKey

  override def coercion = Seq(WdlIntegerType)

  override protected def validateValue = {
    case wdlValue if WdlIntegerType.coerceRawValue(wdlValue).isSuccess =>
      WdlIntegerType.coerceRawValue(wdlValue).get match {
        case WdlInteger(value) => if (value.toInt <= 0) wrongAmountMsg.invalidNel else value.toInt.validNel
      }
  }

  override def validateExpression = {
    case wdlValue if WdlIntegerType.coerceRawValue(wdlValue).isSuccess => true
  }

  override protected def failureMessage = missingMessage
}
