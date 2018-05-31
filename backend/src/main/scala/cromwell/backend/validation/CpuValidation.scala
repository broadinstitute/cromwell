package cromwell.backend.validation

import cats.data.NonEmptyList
import cats.syntax.either._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineV
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
  lazy val instance: RuntimeAttributesValidation[Int Refined Positive] = new CpuValidation(CpuKey)
  lazy val optional: OptionalRuntimeAttributesValidation[Int Refined Positive] = instance.optional
  lazy val instanceMin: RuntimeAttributesValidation[Int Refined Positive] = new CpuValidation(CpuMinKey)
  lazy val optionalMin: OptionalRuntimeAttributesValidation[Int Refined Positive] = instanceMin.optional
  lazy val instanceMax: RuntimeAttributesValidation[Int Refined Positive] = new CpuValidation(CpuMaxKey)
  lazy val optionalMax: OptionalRuntimeAttributesValidation[Int Refined Positive] = instanceMax.optional

  lazy val defaultMin: WomValue = WomInteger(1)
  def configDefaultWomValue(config: Option[Config]): Option[WomValue] = instance.configDefaultWomValue(config)
}

class CpuValidation(attributeName: String) extends PositiveIntRuntimeAttributesValidation(attributeName) {
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Int Refined Positive]] = {
    case womValue if WomIntegerType.coerceRawValue(womValue).isSuccess =>
      WomIntegerType.coerceRawValue(womValue).get match {
        case WomInteger(value) =>
          refineV[Positive](value.toInt)
            .leftMap(_ => NonEmptyList.one(s"Expecting $key runtime attribute value greater than 0"))
            .toValidated
      }
  }

  override protected def missingValueMessage: String = s"Expecting $key runtime attribute to be an Integer"
}
