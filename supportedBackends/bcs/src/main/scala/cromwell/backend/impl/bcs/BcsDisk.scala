package cromwell.backend.impl.bcs

import scala.util.{Try, Success, Failure}

trait BcsDisk {
  val diskType: String
  val sizeInGB: Int
}

final case class BcsSystemDisk(diskType: String, sizeInGB: Int) extends BcsDisk
final case class BcsDataDisk(diskType: String, sizeInGB: Int, mountPoint: String) extends BcsDisk

object BcsDisk{
  val systemDiskPattern = s"""(\\S+)\\s+(\\d+)""".r
  val dataDiskPattern = s"""(\\S+)\\s+(\\d+)\\s+(\\S+)""".r

  def parse(s: String): Try[BcsDisk] = {
    s match {
      case systemDiskPattern(diskType, sizeInGB) => Success(BcsSystemDisk(diskType, sizeInGB.toInt))
      case dataDiskPattern(diskType, sizeInGB, mountPoint) => Success(BcsDataDisk(diskType, sizeInGB.toInt, mountPoint))
      case _ => Failure(new IllegalArgumentException("disk must be 'cloud 40' or 'cloud 200 /home/input/'"))
    }
  }
}