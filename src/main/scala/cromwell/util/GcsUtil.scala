package cromwell.util

import java.nio.file.Path

object GcsUtil {

  /** Read an entire file into a string, closing the underlying stream. */
  def slurp(file: GoogleCloudStoragePath, clientSecretsFile: Path): String = {
    new String(new GcsConnector("cromwell", clientSecretsFile).downloadObject(file), "UTF-8")
  }
}
