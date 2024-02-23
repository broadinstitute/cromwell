package cromwell.backend.google.batch.util

import cats.data.NonEmptyList
import cats.syntax.either._
import cats.syntax.validated._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.validation.{
  OptionalRuntimeAttributesValidation,
  PositiveIntRuntimeAttributesValidation,
  RuntimeAttributesValidation
}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineV
import wom.RuntimeAttributesKeys.GpuKey
import wom.types.WomIntegerType
import wom.values.{WomInteger, WomValue}

object GpuValidation {
  lazy val instance: RuntimeAttributesValidation[Int Refined Positive] = new GpuValidation(GpuKey)
  lazy val optional: OptionalRuntimeAttributesValidation[Int Refined Positive] = instance.optional

  lazy val defaultMin: WomValue = WomInteger(0)
  def configDefaultWomValue(config: Option[Config]): Option[WomValue] = instance.configDefaultWomValue(config)
}

class GpuValidation(attributeName: String) extends PositiveIntRuntimeAttributesValidation(attributeName) {
  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Int Refined Positive]] = {
    case womValue if WomIntegerType.coerceRawValue(womValue).isSuccess =>
      WomIntegerType.coerceRawValue(womValue).get match {
        case WomInteger(value) =>
          refineV[Positive](value.toInt)
            .leftMap(_ => NonEmptyList.one(s"Expecting $key runtime attribute value greater than 0"))
            .toValidated
      }
    case other =>
      s"Invalid gpu count. Expected positive Int but got ${other.womType.friendlyName} ${other.toWomString}".invalidNel
  }
}
