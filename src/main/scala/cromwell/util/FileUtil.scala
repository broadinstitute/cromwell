package cromwell.util

import java.io.{BufferedWriter, Writer, FileWriter, File}
import java.nio.file.Path

object FileUtil {

  /** Build a `Writer` for the specified `File`, return the tuple of both. */
  def fileAndWriter(file: Path): (Path, Writer) = {
    val writer = new BufferedWriter(new FileWriter(file.toFile))
    (file, writer)
  }

  /** Build a temp file with the specified base name and an associated writer,
    * return the tuple of both. */
  def tempFileAndWriter(baseName: String, extension: String = ".tmp"): (Path, Writer) = {
    val file = File.createTempFile(baseName, extension).toPath
    fileAndWriter(file)
  }

  implicit class FlushingAndClosingWriter(writer: Writer) {
    /** Convenience method to flush and close in one shot. */
    def flushAndClose() = {
      writer.flush()
      writer.close()
    }
  }
}
