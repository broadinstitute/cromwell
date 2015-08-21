package cromwell.util

import java.io.{BufferedWriter, File, FileWriter, Writer}
import java.nio.file.Path

import scala.util.{Failure, Success, Try}

object FileUtil {
  /** Build a temp file with the specified base name and an associated writer,
    * return the tuple of both. */
  def tempFileAndWriter(baseName: String, directory: File = null): (Path, Writer) = {
    File.createTempFile(baseName, ".tmp", directory).toPath.fileAndWriter
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

  implicit class EnhancedFile(val file: File) extends AnyVal {
    /** Read an entire file into a string, closing the underlying stream. */
    def slurp: String = {
      val source = io.Source.fromFile(file)
      try source.mkString finally source.close()
    }
  }

  implicit class EnhancedPath(val path: Path) extends AnyVal {
    def slurp = path.toFile.slurp

    def fileAndWriter: (Path, Writer) = {
      val writer = new BufferedWriter(new FileWriter(path.toFile))
      (path, writer)
    }
  }
}
