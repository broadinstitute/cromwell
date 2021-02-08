package cromwell.backend.google.pipelines.common.errors

import scala.util.control.NoStackTrace

class InvalidGcsPathsInManifestFileException(paths: List[String]) extends Exception with NoStackTrace {
  override def getMessage: String = s"Some of the paths in manifest file are not valid GCS paths: \n${paths.mkString("\n")}"
}
