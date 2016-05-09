package cromwell.util

import java.io.{File, FileInputStream, Writer}

import better.files._
import org.apache.commons.codec.digest.DigestUtils
import wdl4s.values.Hashable

import scala.util.{Failure, Success, Try}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
object FileUtil {
  def swapExt(filePath: String, oldExt: String, newExt: String): String = {
    filePath.stripSuffix(oldExt) + newExt
  }

  def parseTsv(tsv: String): Try[Array[Array[String]]] = {
    val table = tsv.split("\n").map(_.split("\t"))
    table.map(_.size).toSet match {
      case s if s.size > 1 => Failure(new UnsupportedOperationException("TSV is not uniform"))
      case _ => Success(table)
    }
  }

  implicit class FlushingAndClosingWriter(writer: Writer) {
    /** Convenience method to flush and close in one shot. */
    def flushAndClose() = {
      writer.flush()
      writer.close()
    }
  }

  implicit class EnhancedFile(val file: File) extends AnyVal with Hashable {
    /** Read an entire file into a string, closing the underlying stream. */
    def slurp: String = {
      // TODO: deprecate slurp, and java.io.File in general?
      file.toPath.contentAsString
    }

    def md5Sum: String = {
      val fis = new FileInputStream(file)
      try {
        DigestUtils.md5Hex(fis)
      } finally fis.close()
    }
  }
}