package cromwell.backend.google.batch.io

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import cromwell.backend.DiskPatterns._
import cromwell.core.path.{DefaultPathBuilder, Path}
import wom.values._

import scala.util.Try

object GcpBatchAttachedDisk {
  def parse(s: String): Try[GcpBatchAttachedDisk] = {

    def sizeGbValidation(sizeGbString: String): ErrorOr[Int] = validateLong(sizeGbString).map(_.toInt)

    def diskTypeValidation(diskTypeString: String): ErrorOr[DiskType] = validateDiskType(diskTypeString)

    val validation: ErrorOr[GcpBatchAttachedDisk] = s match {
      case WorkingDiskPattern(sizeGb, diskType) =>
        (validateDiskType(diskType), sizeGbValidation(sizeGb)) mapN {
          GcpBatchWorkingDisk.apply
        }
      case MountedDiskPattern(mountPoint, sizeGb, diskType) =>
        (sizeGbValidation(sizeGb), diskTypeValidation(diskType)) mapN { (s, dt) =>
          PipelinesApiEmptyMountedDisk(dt, s, DefaultPathBuilder.get(mountPoint))
        }
      case _ =>
        s"Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE' but got: '$s'".invalidNel
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

  private def validateDiskType(diskTypeName: String): ErrorOr[DiskType] =
    DiskType.values().find(_.diskTypeName == diskTypeName) match {
      case Some(diskType) => diskType.validNel
      case None =>
        val diskTypeNames = DiskType.values.map(_.diskTypeName).mkString(", ")
        s"Disk TYPE $diskTypeName should be one of $diskTypeNames".invalidNel
    }

  private def validateLong(value: String): ErrorOr[Long] =
    try
      value.toLong.validNel
    catch {
      case _: IllegalArgumentException => s"$value not convertible to a Long".invalidNel
    }

}

trait GcpBatchAttachedDisk {
  def name: String
  def diskType: DiskType
  def sizeGb: Int
  def mountPoint: Path
}

case class PipelinesApiEmptyMountedDisk(diskType: DiskType, sizeGb: Int, mountPoint: Path)
    extends GcpBatchAttachedDisk {
  val name = s"d-${mountPoint.pathAsString.md5Sum}"

  override def toString: String = s"$mountPoint $sizeGb ${diskType.diskTypeName}"
}

object GcpBatchWorkingDisk {
  val MountPoint: Path = DefaultPathBuilder.get("/mnt/disks/cromwell_root")
  val Name = "local-disk"
  val Default = GcpBatchWorkingDisk(DiskType.SSD, 10)
}

case class GcpBatchWorkingDisk(diskType: DiskType, sizeGb: Int) extends GcpBatchAttachedDisk {
  val mountPoint: Path = GcpBatchWorkingDisk.MountPoint
  val name: String = GcpBatchWorkingDisk.Name

  override def toString: String = s"$name $sizeGb ${diskType.diskTypeName}"
}

case class GcpBatchReferenceFilesDisk(image: String, sizeGb: Int) extends GcpBatchAttachedDisk {
  val mountPoint: Path = DefaultPathBuilder.get(s"/mnt/${image.md5Sum}")
  val name: String = s"d-${mountPoint.pathAsString.md5Sum}"
  val diskType: DiskType = DiskType.HDD
}
