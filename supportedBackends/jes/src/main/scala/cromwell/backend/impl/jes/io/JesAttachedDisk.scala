package cromwell.backend.impl.jes.io

import java.nio.file.{Path, Paths}

import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._
import com.google.api.services.genomics.model.Disk
import lenthall.validation.ErrorOr._
import lenthall.exception.MessageAggregation
import wdl4s.values._

import scala.util.Try


object JesAttachedDisk {
  val Identifier = "[a-zA-Z0-9-_]+"
  val Directory = """/[^\s]+"""
  val Integer = "[1-9][0-9]*"
  val WorkingDiskPattern = s"""${JesWorkingDisk.Name}\\s+($Integer)\\s+($Identifier)""".r
  val MountedDiskPattern = s"""($Directory)\\s+($Integer)\\s+($Identifier)""".r

  def parse(s: String): Try[JesAttachedDisk] = {
    val validation: ErrorOr[JesAttachedDisk] = s match {
      case WorkingDiskPattern(sizeGb, diskType) =>
        (validateLong(sizeGb) |@| validateDiskType(diskType)) map { (s, dt) =>
          JesWorkingDisk(dt, s.toInt)
        }
      case MountedDiskPattern(mountPoint, sizeGb, diskType) =>
        (validateLong(sizeGb) |@| validateDiskType(diskType)) map { (s, dt) =>
          JesEmptyMountedDisk(dt, s.toInt, Paths.get(mountPoint))
        }
      case _ => s"Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE'".invalidNel
    }

    Try(validation match {
      case Valid(localDisk) => localDisk
      case Invalid(nels) =>
        throw new UnsupportedOperationException with MessageAggregation {
          val exceptionContext = ""
          val errorMessages = nels.toList
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
}

trait JesAttachedDisk {
  def name: String
  def diskType: DiskType
  def sizeGb: Int
  def mountPoint: Path
  def toGoogleDisk: Disk = {
    new Disk().setName(name)
      .setType(diskType.googleTypeName)
      .setAutoDelete(true)
      .setSizeGb(sizeGb)
      .setMountPoint(mountPoint.toAbsolutePath.toString)
  }
}

case class JesEmptyMountedDisk(diskType: DiskType, sizeGb: Int, mountPoint: Path) extends JesAttachedDisk {
  val name = s"d-${mountPoint.toString.md5Sum}"
  override def toString: String = s"$mountPoint $sizeGb ${diskType.diskTypeName}"
}

object JesWorkingDisk {
  val MountPoint = Paths.get("/cromwell_root")
  val Name = "local-disk"
}

case class JesWorkingDisk(diskType: DiskType, sizeGb: Int) extends JesAttachedDisk {
  val mountPoint = JesWorkingDisk.MountPoint
  val name = JesWorkingDisk.Name
  override def toString: String = s"$name $sizeGb ${diskType.diskTypeName}"
}
