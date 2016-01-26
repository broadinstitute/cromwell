package cromwell.engine.backend.runtimeattributes

import scala.util.Try
import scalaz._
import scalaz.Scalaz._
import com.google.api.services.genomics.model.Disk

object LocalDisk {
  def parse(s: String): Try[LocalDisk] = {
    val validation = s.split("\\s+") match {
      case Array(name, DiskType.LOCAL.diskTypeName) =>
        LocalDisk(name, DiskType.LOCAL).successNel[String]
      case Array(name, sizeGb, diskType) if diskType != DiskType.LOCAL.diskTypeName =>
        (validateLong(sizeGb) |@| validateDiskType(diskType)) { (s, dt) => LocalDisk(name, dt, s) }
      case _ => s"'$s' should be in form 'NAME SIZE TYPE', with SIZE blank for LOCAL, otherwise SIZE in GB".failureNel
    }

    Try(validation match {
      case Success(localDisk) => localDisk
      case Failure(nels) => throw new UnsupportedOperationException(nels.list.mkString("\n"))
    })
  }

  private def validateDiskType(diskTypeName: String): ValidationNel[String, DiskType] = {
    DiskType.values().find(_.diskTypeName == diskTypeName) match {
      case Some(diskType) => diskType.successNel[String]
      case None =>
        val diskTypeNames = DiskType.values.map(_.diskTypeName).mkString(", ")
        s"Disk TYPE $diskTypeName should be one of $diskTypeNames".failureNel
    }
  }

  private def validateLong(value: String): ValidationNel[String, Long] = {
    try {
      value.toLong.successNel
    } catch {
      case _: IllegalArgumentException => s"$value not convertible to a Long".failureNel[Long]
    }
  }
}

case class LocalDisk(name: String, diskType: DiskType, sizeGb: Long = 10L) {
  def toDisk: Disk = {
    val disk = new Disk().setName(name).setType(diskType.googleTypeName).setAutoDelete(true)
    // Even though GCE ignores the value for local disks, JES requires we set the disk size anyway.
    disk.setSizeGb(long2Long(sizeGb))
    disk
  }
  override def toString: String = s"$name $sizeGb $diskType"
}

