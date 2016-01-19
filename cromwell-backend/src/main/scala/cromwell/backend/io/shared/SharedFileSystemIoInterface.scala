package cromwell.backend.io.shared

import java.nio.file.{Files, Paths}

import cromwell.backend.io.IoInterface

import scala.language.postfixOps
import scala.util.Try

object PathUtil {

  implicit class UriString(val str: String) extends AnyVal {
    def isGcsUrl: Boolean = str.startsWith("gs://")
    def isUriWithProtocol: Boolean = "^[a-z]+://".r.findFirstIn(str).nonEmpty
  }

}

object SharedFileSystemIoInterface {
  lazy val instance = new SharedFileSystemIoInterface()
}

class SharedFileSystemIoInterface private() extends IoInterface {
  import PathUtil._
  import better.files._

  override def readFile(path: String): String = Paths.get(path).contentAsString

  override def writeFile(path: String, content: String): Unit = Paths.get(path).write(content)

  override def listContents(path: String): Iterable[String] = Paths.get(path).list map {
    _.path.toAbsolutePath.toString
  } toIterable

  override def exists(path: String): Boolean = Paths.get(path).exists

  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = {
    val file = Files.createTempFile(Paths.get(path), prefix, suffix)
    file.write(content)
    file.fullPath
  }

  override def glob(path: String, pattern: String): Seq[String] = {
    Paths.get(path).glob(s"**/$pattern") map { _.path.fullPath } toSeq
  }

  def isValidPath(path: String) = !path.isUriWithProtocol && Try(Paths.get(path)).isSuccess

  override def copy(from: String, to: String): Unit = {
    val src = new File(Paths.get(from))
    val dst = new File(Paths.get(to))
    src.copyTo(dst)
  }

  override def hash(path: String) = new File(Paths.get(path)).md5
}
