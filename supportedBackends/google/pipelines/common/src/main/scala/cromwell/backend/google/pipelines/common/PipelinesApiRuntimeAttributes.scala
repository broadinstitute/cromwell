package cromwell.backend.google.pipelines.common

import cats.data.Validated._
import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import cromwell.backend.google.pipelines.common.io.{PipelinesApiAttachedDisk, PipelinesApiWorkingDisk}
import cromwell.backend.validation._
import wom.types._
import wom.values._

object DisksValidation extends RuntimeAttributesValidation[Seq[PipelinesApiAttachedDisk]] {
  override def key: String = "disks"

  override def coercion: Iterable[WomType] = Set(WomStringType, WomArrayType(WomStringType))

  override protected def validateValue: PartialFunction[WomValue, ErrorOr[Seq[PipelinesApiAttachedDisk]]] = {
    case WomString(value) => validateLocalDisks(value.split(",\\s*").toSeq)
    case WomArray(womType, values) if womType.memberType == WomStringType =>
      validateLocalDisks(values.map(_.valueString))
  }

  private def validateLocalDisks(disks: Seq[String]): ErrorOr[Seq[PipelinesApiAttachedDisk]] = {
    val diskNels: ErrorOr[Seq[PipelinesApiAttachedDisk]] =
      disks.toList.traverse[ErrorOr, PipelinesApiAttachedDisk](validateLocalDisk)
    val defaulted: ErrorOr[Seq[PipelinesApiAttachedDisk]] = addDefault(diskNels)
    defaulted
  }

  private def validateLocalDisk(disk: String): ErrorOr[PipelinesApiAttachedDisk] =
    PipelinesApiAttachedDisk.parse(disk) match {
      case scala.util.Success(attachedDisk) => attachedDisk.validNel
      case scala.util.Failure(ex) => ex.getMessage.invalidNel
    }

  private def addDefault(disksNel: ErrorOr[Seq[PipelinesApiAttachedDisk]]): ErrorOr[Seq[PipelinesApiAttachedDisk]] =
    disksNel map {
      case disks if disks.exists(_.name == PipelinesApiWorkingDisk.Name) => disks
      case disks => disks :+ PipelinesApiWorkingDisk.Default
    }

  override protected def missingValueMessage: String =
    s"Expecting $key runtime attribute to be a comma separated String or Array[String]"
}
