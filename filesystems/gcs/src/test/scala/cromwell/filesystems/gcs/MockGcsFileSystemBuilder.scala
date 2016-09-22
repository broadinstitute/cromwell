package cromwell.filesystems.gcs

import scala.util.Failure

object MockGcsFileSystemBuilder {
  val mockGcsFileSystem = new GcsFileSystemProvider(
    Failure(new Exception("No Storage object available")),
    scala.concurrent.ExecutionContext.global).defaultFileSystem
}
