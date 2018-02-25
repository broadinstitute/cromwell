package cromwell.backend.impl.bcs

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import common.exception.MessageAggregation
import common.validation.ErrorOr._
import cromwell.core.path.{DefaultPathBuilder, Path, PathBuilder, PathFactory}

import scala.util.Try
import scala.util.matching.Regex

object BcsMount {
  var pathBuilders: List[PathBuilder] = List()

  val ossPrefix = """oss://[^\s]+"""
  val localPath = """/[^\s]+"""
  val writeSupport = """true|false"""

  val inputMountPattern: Regex = s"""($ossPrefix)\\s+($localPath)\\s+($writeSupport)""".r
  val outputMountPattern: Regex = s"""($localPath)\\s+($ossPrefix)\\s+($writeSupport)""".r

  def parse(s: String): Try[BcsMount] = {
    val validation: ErrorOr[BcsMount] = s match {
      case inputMountPattern(oss, local, writeSupport) =>
        (validateOss(oss), validateLocal(oss, local), validateBoolean(writeSupport)) mapN { (src, dest, ws) => new BcsInputMount(src, dest, ws)}
      case outputMountPattern(local, oss, writeSupport) =>
        (validateLocal(oss, local), validateOss(oss), validateBoolean(writeSupport)) mapN { (src, dest, ws) => new BcsOutputMount(src, dest, ws)}
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

  private def validateOss(value: String): ErrorOr[Path] = {
    PathFactory.buildPath(value, pathBuilders).validNel
  }
  private def validateLocal(oss: String, local: String): ErrorOr[Path] = {
    if (oss.endsWith("/") == local.endsWith("/")) {
      DefaultPathBuilder.get(local).validNel
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
  var src: Path
  var dest: Path
  var writeSupport: Boolean
}

final case class BcsInputMount(var src: Path, var dest: Path, var writeSupport: Boolean) extends BcsMount
final case class BcsOutputMount(var src: Path, var dest: Path, var writeSupport: Boolean) extends BcsMount