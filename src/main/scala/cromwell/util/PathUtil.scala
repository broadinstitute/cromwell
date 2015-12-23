package cromwell.util

import java.nio.file.{Path, Paths}

import cromwell.engine.io.gcs.GcsFileSystem


object PathUtil {

  implicit class UriString(val str: String) extends AnyVal {
    def isGcsUrl: Boolean = str.startsWith("gs://")
    def isUriWithProtocol: Boolean = "^[a-z]+://".r.findFirstIn(str).nonEmpty

    def toPath(implicit gcsFileSystem: Option[GcsFileSystem] = None): Path = {
      str match {
        case path if path.isGcsUrl && gcsFileSystem.isDefined => gcsFileSystem.get.getPath(str)
        case path if !path.isUriWithProtocol => Paths.get(path)
        case path => throw new Throwable(s"Unable to parse $path")
      }
    }
  }

}
