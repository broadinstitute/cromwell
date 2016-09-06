package cromwell.util

import java.nio.file.Path

import better.files._
import wdl4s.values.Hashable

import scala.util.{Failure, Success, Try}

object FileUtil {
  def swapExt(filePath: String, oldExt: String, newExt: String): String = {
    filePath.stripSuffix(oldExt) + newExt
  }

  def parseTsv(tsv: String): Try[Array[Array[String]]] = {
    val table = tsv.split("\n").map(_.split("\t"))
    table.map(_.length).toSet match {
      case s if s.size > 1 => Failure(new UnsupportedOperationException("TSV is not uniform"))
      case _ => Success(table)
    }
  }

  implicit class EnhancedFile(val file: Path) extends AnyVal with Hashable {
    def md5Sum: String = File(file).md5.toLowerCase // toLowerCase for backwards compatibility
  }
}
