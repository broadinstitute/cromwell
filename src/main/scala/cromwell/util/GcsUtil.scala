package cromwell.util

object GcsUtil {

  /** Read an entire file into a string, closing the underlying stream. */
  def slurp(file: GoogleCloudStoragePath, clientSecretsFile: String): String = {
    new String(new GcsConnector("cromwell", clientSecretsFile).downloadObject(file), "UTF-8")
  }
}
