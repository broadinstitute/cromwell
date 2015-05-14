package cromwell.util

import java.io.{BufferedWriter, Writer, FileWriter, File}

object FileUtil {

  /** Build a `Writer` for the specified `File`, return the tuple of both. */
  def fileAndWriter(file: File): (File, Writer) = {
    val writer = new BufferedWriter(new FileWriter(file))
    (file, writer)
  }

  /** Build a temp file with the specified base name and an associated writer,
    * return the tuple of both. */
  def tempFileAndWriter(baseName: String, extension: String = ".tmp"): (File, Writer) = {
    val file = File.createTempFile(baseName, extension)
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
