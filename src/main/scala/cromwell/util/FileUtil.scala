package cromwell.util

import java.io.{BufferedWriter, File, FileWriter, Writer}
import java.nio.file.Path

object FileUtil {

  /** Build a `Writer` for the specified `File`, return the tuple of both. */
  def fileAndWriter(file: Path): (Path, Writer) = {
    val writer = new BufferedWriter(new FileWriter(file.toFile))
    (file, writer)
  }

  /** Build a temp file with the specified base name and an associated writer,
    * return the tuple of both. */
  def tempFileAndWriter(baseName: String, directory: File = null): (Path, Writer) = {
    val file = File.createTempFile(baseName, ".tmp", directory).toPath
    fileAndWriter(file)
  }

  implicit class FlushingAndClosingWriter(writer: Writer) {
    /** Convenience method to flush and close in one shot. */
    def flushAndClose() = {
      writer.flush()
      writer.close()
    }
  }

  /** Read an entire file into a string, closing the underlying stream. */
  def slurp(path: Path): String = slurp(path.toFile)

  /** Read an entire file into a string, closing the underlying stream. */
  def slurp(file: File): String = {
    val source = io.Source.fromFile(file)
    try source.mkString finally source.close()
  }
}
