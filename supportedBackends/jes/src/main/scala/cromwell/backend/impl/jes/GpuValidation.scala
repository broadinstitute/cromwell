package cromwell.backend.impl.jes

import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.validation.{IntRuntimeAttributesValidation, OptionalRuntimeAttributesValidation, RuntimeAttributesValidation}
import wom.RuntimeAttributesKeys.{GpuKey, GpuMaxKey, GpuMinKey}
import wom.types.WomIntegerType
import wom.values.{WomInteger, WomValue}

object GpuValidation {
  lazy val instance: RuntimeAttributesValidation[Int] = new GpuValidation(GpuKey)
  lazy val optional: OptionalRuntimeAttributesValidation[Int] = instance.optional
  lazy val instanceMin: RuntimeAttributesValidation[Int] = new GpuValidation(GpuMinKey)
  lazy val optionalMin: OptionalRuntimeAttributesValidation[Int] = instanceMin.optional
  lazy val instanceMax: RuntimeAttributesValidation[Int] = new GpuValidation(GpuMaxKey)
  lazy val optionalMax: OptionalRuntimeAttributesValidation[Int] = instanceMax.optional

  lazy val defaultMin: WomValue = WomInteger(0)
  def configDefaultWomValue(config: Option[Config]): Option[WomValue] = instance.configDefaultWomValue(config)
}

class GpuValidation(attributeName: String) extends IntRuntimeAttributesValidation(attributeName) {
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Int]] = {
    case womValue if WomIntegerType.coerceRawValue(womValue).isSuccess =>
      WomIntegerType.coerceRawValue(womValue).get match {
        case WomInteger(value) =>
          if (value.toInt < 0)
            s"Expecting $key runtime attribute value greater or equal than 0".invalidNel
          else
            value.toInt.validNel
      }
  }
}
