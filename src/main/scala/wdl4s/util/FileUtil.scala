package wdl4s.util

import java.io.File
import better.files._
import scala.util.{Failure, Success, Try}

object FileUtil {
  // FIXME: Keep?
  def parseTsv(tsv: String): Try[Array[Array[String]]] = {
    val table = tsv.split("\n").map(_.split("\t"))
    table.map(_.size).toSet match {
      case s if s.size > 1 => Failure(new UnsupportedOperationException("TSV is not uniform"))
      case _ => Success(table)
    }
  }

  implicit class EnhancedFile(val file: File) extends AnyVal {
    /** Read an entire file into a string, closing the underlying stream. */
    def slurp: String = {
      // TODO: deprecate slurp, and java.io.File in general?
      file.toPath.contentAsString
    }
  }
}