package cromwell.backend.google.batch.util

import cats.implicits.catsSyntaxValidatedId
import common.validation.ErrorOr.ErrorOr
import cromwell.backend.google.batch.models.{GcpBatchRuntimeAttributes, MachineType}
import cromwell.backend.validation.{OptionalRuntimeAttributesValidation, RuntimeAttributesValidation}
import wom.types.{WomStringType, WomType}
import wom.values.{WomString, WomValue}

object MachineTypeValidation {
  lazy val instance: RuntimeAttributesValidation[MachineType] = new MachineTypeValidation
  lazy val optional: OptionalRuntimeAttributesValidation[MachineType] = instance.optional
}

class MachineTypeValidation extends RuntimeAttributesValidation[MachineType] {
  override def key = GcpBatchRuntimeAttributes.MachineTypeKey

  override def coercion: Iterable[WomType] = Set(WomStringType)
  override def validateValue: PartialFunction[WomValue, ErrorOr[MachineType]] = {
    case WomString(s) => MachineType(s).validNel
    case other =>
      s"Invalid '$key': String value required but got ${other.womType.friendlyName}.".invalidNel
  }
}
