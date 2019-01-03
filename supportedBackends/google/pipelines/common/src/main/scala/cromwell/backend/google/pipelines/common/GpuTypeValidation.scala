package cromwell.backend.google.pipelines.common

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.google.pipelines.common.GpuResource.GpuType
import cromwell.backend.google.pipelines.common.GpuResource.GpuType.GpuType
import cromwell.backend.google.pipelines.common.GpuTypeValidation._
import cromwell.backend.validation.{OptionalRuntimeAttributesValidation, RuntimeAttributesValidation}
import wom.RuntimeAttributesKeys
import wom.types.{WomStringType, WomType}
import wom.values.{WomString, WomValue}

import scala.util.{Failure, Success, Try}

object GpuTypeValidation {
  lazy val instance: RuntimeAttributesValidation[GpuType] = new GpuTypeValidation
  lazy val optional: OptionalRuntimeAttributesValidation[GpuType] = instance.optional
  private val SupportedTypesMessage = s"Supported types are ${GpuResource.GpuType.NVIDIATeslaK80.toString}, ${GpuResource.GpuType.NVIDIATeslaP100.toString}"
}

class GpuTypeValidation extends RuntimeAttributesValidation[GpuType] {
  override def key = RuntimeAttributesKeys.GpuTypeKey

  override def coercion: Traversable[WomType] = Set(WomStringType)
  override def validateValue: PartialFunction[WomValue, ErrorOr[GpuType]] = {
    case WomString(s) => Try(GpuType.withName(s.toString)) match {
      case Success(tpe: GpuType) => tpe.validNel
      case Failure(_: NoSuchElementException) => s"$s is not a supported GPU type. $SupportedTypesMessage".invalidNel
      case Failure(e) => s"Could not parse $s as a valid GPU type. $SupportedTypesMessage. Error: ${e.getMessage}".invalidNel
    }
  }
}
