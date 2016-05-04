package cromwell.engine.backend.jes

import java.nio.file.{Path, Paths}

import com.google.api.services.genomics.model.Disk
import cromwell.core.ErrorOr
import cromwell.engine.backend.runtimeattributes.DiskType
import wdl4s.ThrowableWithErrors
import wdl4s.values._

import scala.util.Try
import scalaz.Scalaz._
import scalaz._

object JesAttachedDisk {
  val Identifier = "[a-zA-Z0-9-_]+"
  val Directory = """/[^\s]+"""
  val Integer = "[1-9][0-9]*"
  val WorkingDiskPattern = s"""${JesWorkingDisk.Name}\\s+($Integer)\\s+($Identifier)""".r
  val MountedDiskPattern = s"""($Directory)\\s+($Integer)\\s+($Identifier)""".r

  def parse(s: String): Try[JesAttachedDisk] = {
    val validation = s match {
      case WorkingDiskPattern(sizeGb, diskType) =>
        (validateLong(sizeGb) |@| validateDiskType(diskType)) { (s, dt) =>
          JesWorkingDisk(dt, s.toInt)
        }
      case MountedDiskPattern(mountPoint, sizeGb, diskType) =>
        (validateLong(sizeGb) |@| validateDiskType(diskType)) { (s, dt) =>
          JesEmptyMountedDisk(dt, s.toInt, Paths.get(mountPoint))
        }
      case _ => s"Disk strings should be of the format 'local-disk SIZE TYPE' or '/mount/point SIZE TYPE'".failureNel
    }

    Try(validation match {
      case Success(localDisk) => localDisk
      case Failure(nels) =>
        throw new UnsupportedOperationException with ThrowableWithErrors {
          val message = ""
          val errors = nels
        }
    })
  }

  private def validateDiskType(diskTypeName: String): ErrorOr[DiskType] = {
    DiskType.values().find(_.diskTypeName == diskTypeName) match {
      case Some(diskType) => diskType.successNel[String]
      case None =>
        val diskTypeNames = DiskType.values.map(_.diskTypeName).mkString(", ")
        s"Disk TYPE $diskTypeName should be one of $diskTypeNames".failureNel
    }
  }

  private def validateLong(value: String): ErrorOr[Long] = {
    try {
      value.toLong.successNel
    } catch {
      case _: IllegalArgumentException => s"$value not convertible to a Long".failureNel[Long]
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
  val MountPoint = "/cromwell_root"
  val Name = "local-disk"
}

case class JesWorkingDisk(diskType: DiskType, sizeGb: Int) extends JesAttachedDisk {
  val mountPoint = Paths.get(JesWorkingDisk.MountPoint)
  val name = JesWorkingDisk.Name
  override def toString: String = s"$name $sizeGb ${diskType.diskTypeName}"
}
