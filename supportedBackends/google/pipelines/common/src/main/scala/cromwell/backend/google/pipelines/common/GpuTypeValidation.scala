package cromwell.backend.google.pipelines.common

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.google.pipelines.common.GpuResource.GpuType
import cromwell.backend.validation.{OptionalRuntimeAttributesValidation, RuntimeAttributesValidation}
import wom.RuntimeAttributesKeys
import wom.types.{WomStringType, WomType}
import wom.values.{WomString, WomValue}

object GpuTypeValidation {
  lazy val instance: RuntimeAttributesValidation[GpuType] = new GpuTypeValidation
  lazy val optional: OptionalRuntimeAttributesValidation[GpuType] = instance.optional
}

class GpuTypeValidation extends RuntimeAttributesValidation[GpuType] {
  override def key = RuntimeAttributesKeys.GpuTypeKey

  override def coercion: Traversable[WomType] = Set(WomStringType)
  override def validateValue: PartialFunction[WomValue, ErrorOr[GpuType]] = {
    case WomString(s) => GpuType(s).validNel
    case other => s"Invalid '$key': String value required but got ${other.womType.friendlyName}. See ${GpuType.MoreDetailsURL} for a list of options".invalidNel
  }
}
