package cromwell.backend.impl.bcs

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import com.aliyuncs.batchcompute.pojo.v20151111.MountEntry
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import cromwell.backend.impl.bcs.BcsMount.PathType
import cromwell.core.path.{Path, PathBuilder, PathFactory}

import scala.util.{Success, Try}
import scala.util.matching.Regex

object BcsMount {
  type PathType = Either[Path, String]

  def toString(p: PathType): String = {
    p match {
      case Left(p) =>
        p.pathAsString
      case Right(s) =>
        return s
    }
  }

  val supportFileSystemTypes = List("oss", "nas", "smb", "lustre").mkString("|")

  var pathBuilders: List[PathBuilder] = List()

  val remotePrefix = s"""(?:$supportFileSystemTypes)""" + """://[^\s]+"""
  val localPath = """/[^\s]+"""
  val writeSupport = """true|false"""

  val inputMountPattern: Regex = s"""($remotePrefix)\\s+($localPath)\\s+($writeSupport)""".r
  val outputMountPattern: Regex = s"""($localPath)\\s+($remotePrefix)\\s+($writeSupport)""".r

  def parse(s: String): Try[BcsMount] = {
    val validation: ErrorOr[BcsMount] = s match {
      case inputMountPattern(remote, local, writeSupport) =>
        (validateRemote(remote), validateLocal(remote, local), validateBoolean(writeSupport)) mapN { (src, dest, ws) => new BcsInputMount(src, dest, ws)}
      case outputMountPattern(local, oss, writeSupport) =>
        (validateLocal(oss, local), validateRemote(oss), validateBoolean(writeSupport)) mapN { (src, dest, ws) => new BcsOutputMount(src, dest, ws)}
      case _ => s"Mount strings should be of the format 'oss://my-bucket/inputs/ /home/inputs/ true' or '/home/outputs/ oss://my-bucket/outputs/ false'".invalidNel
    }

    Try(validation match {
      case Valid(mount) => mount
      case Invalid(nels) =>
        throw new UnsupportedOperationException with MessageAggregation {
          val exceptionContext = ""
          val errorMessages: List[String] = nels.toList
        }
    })
  }

  private def validateRemote(value: String): ErrorOr[PathType] = {
    Try(PathFactory.buildPath(value, pathBuilders)) match {
      case Success(p) =>
        Left(p).validNel
      case _ =>
        Right(value).validNel
    }
  }
  private def validateLocal(remote: String, local: String): ErrorOr[PathType] = {
    if (remote.endsWith("/") == local.endsWith("/")) {
      Try(PathFactory.buildPath(local, pathBuilders)) match {
        case Success(p) =>
          Left(p).validNel
        case _=>
          Right(local).validNel
      }
    } else {
      "oss and local path type not match".invalidNel
    }
  }

  private def validateBoolean(value: String): ErrorOr[Boolean] = {
    try {
      value.toBoolean.validNel
    } catch {
      case _: IllegalArgumentException => s"$value not convertible to a Boolean".invalidNel
    }
  }
}

trait BcsMount {
  import BcsMount._
  var src: PathType
  var dest: PathType
  var writeSupport: Boolean

  def toBcsMountEntry: MountEntry
}

final case class BcsInputMount(var src: PathType, var dest: PathType, var writeSupport: Boolean) extends BcsMount {
  def toBcsMountEntry: MountEntry = {
    var destStr = BcsMount.toString(dest)
    if (BcsMount.toString(src).endsWith("/") && !destStr.endsWith("/")) {
      destStr += "/"
    }

    val entry = new MountEntry
    entry.setSource(BcsMount.toString(src))
    entry.setDestination(destStr)
    entry.setWriteSupport(writeSupport)

    entry
  }

}
final case class BcsOutputMount(var src: PathType, var dest: PathType, var writeSupport: Boolean) extends BcsMount {
  def toBcsMountEntry: MountEntry = {
    var srcStr = BcsMount.toString(src)
    if (BcsMount.toString(dest).endsWith("/") && !srcStr.endsWith("/")) {
      srcStr += "/"
    }


    val entry = new MountEntry
    entry.setSource(srcStr)
    entry.setDestination(BcsMount.toString(dest))
    entry.setWriteSupport(writeSupport)

    entry
  }
}