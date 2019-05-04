package cromwell.backend.google.pipelines.common.io

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import cromwell.backend.DiskPatterns._
import cromwell.core.path.{DefaultPathBuilder, Path}
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.values._

import scala.util.Try

object PipelinesApiAttachedDisk {
  def parse(s: String): Try[PipelinesApiAttachedDisk] = {

    def sizeGbValidation(sizeGbString: String): ErrorOr[Int] = validateLong(sizeGbString).map(_.toInt)
    def diskTypeValidation(diskTypeString: String): ErrorOr[DiskType] = validateDiskType(diskTypeString)

    val validation: ErrorOr[PipelinesApiAttachedDisk] = s match {
      case WorkingDiskPattern(sizeGb, diskType) => (validateDiskType(diskType), sizeGbValidation(sizeGb)) mapN { PipelinesApiWorkingDisk.apply }
      case MountedDiskPattern(mountPoint, sizeGb, diskType) => (sizeGbValidation(sizeGb), diskTypeValidation(diskType)) mapN { (s, dt) => PipelinesApiEmptyMountedDisk(dt, s, DefaultPathBuilder.get(mountPoint)) }
      case _ => s"Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE' but got: '$s'".invalidNel
    }

    Try(validation match {
      case Valid(localDisk) => localDisk
      case Invalid(nels) =>
        throw new UnsupportedOperationException with MessageAggregation {
          val exceptionContext = ""
          val errorMessages: List[String] = nels.toList
        }
    })
  }

  private def validateDiskType(diskTypeName: String): ErrorOr[DiskType] = {
    DiskType.values().find(_.diskTypeName == diskTypeName) match {
      case Some(diskType) => diskType.validNel
      case None =>
        val diskTypeNames = DiskType.values.map(_.diskTypeName).mkString(", ")
        s"Disk TYPE $diskTypeName should be one of $diskTypeNames".invalidNel
    }
  }

  private def validateLong(value: String): ErrorOr[Long] = {
    try {
      value.toLong.validNel
    } catch {
      case _: IllegalArgumentException => s"$value not convertible to a Long".invalidNel
    }
  }
  
  implicit class EnhancedDisks(val disks: Seq[PipelinesApiAttachedDisk]) extends AnyVal {
    def adjustWorkingDiskWithNewMin(minimum: MemorySize, onAdjustment: => Unit): Seq[PipelinesApiAttachedDisk] = disks map {
      case disk: PipelinesApiWorkingDisk if disk == PipelinesApiWorkingDisk.Default && disk.sizeGb < minimum.to(MemoryUnit.GB).amount.toInt =>
        onAdjustment
        disk.copy(sizeGb = minimum.to(MemoryUnit.GB).amount.toInt)
      case other => other
    }
  }
}

trait PipelinesApiAttachedDisk {
  def name: String
  def diskType: DiskType
  def sizeGb: Int
  def mountPoint: Path
}

case class PipelinesApiEmptyMountedDisk(diskType: DiskType, sizeGb: Int, mountPoint: Path) extends PipelinesApiAttachedDisk {
  val name = s"d-${mountPoint.pathAsString.md5Sum}"
  override def toString: String = s"$mountPoint $sizeGb ${diskType.diskTypeName}"
}

object PipelinesApiWorkingDisk {
  val MountPoint: Path = DefaultPathBuilder.get("/cromwell_root")
  val Name = "local-disk"
  val Default = PipelinesApiWorkingDisk(DiskType.SSD, 10)
}

case class PipelinesApiWorkingDisk(diskType: DiskType, sizeGb: Int) extends PipelinesApiAttachedDisk {
  val mountPoint = PipelinesApiWorkingDisk.MountPoint
  val name = PipelinesApiWorkingDisk.Name
  override def toString: String = s"$name $sizeGb ${diskType.diskTypeName}"
}
