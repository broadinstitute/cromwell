package cromwell.backend.impl.jes

import java.nio.file.{FileSystem, Path}

import cromwell.core.{PathFactory, WorkflowOptions}
import cromwell.filesystems.gcs.GoogleAuthMode.GoogleAuthOptions
import cromwell.filesystems.gcs.{GcsFileSystem, GoogleAuthMode}

import scala.util.Try

object JesImplicits {
  implicit class GoogleAuthWorkflowOptions(val workflowOptions: WorkflowOptions) extends AnyVal {
    def toGoogleAuthOptions: GoogleAuthMode.GoogleAuthOptions = new GoogleAuthOptions {
      override def get(key: String): Try[String] = workflowOptions.get(key)
    }
  }

  object PathBuilder extends PathFactory

  implicit class PathString(val str: String) extends AnyVal {
    def isGcsUrl: Boolean = str.startsWith("gs://")
    def isUriWithProtocol: Boolean = "^[a-z]+://".r.findFirstIn(str).nonEmpty

    def toPath(fss: List[FileSystem]): Path = PathBuilder.buildPath(str, fss)
    def toPath(fs: FileSystem): Path = str.toPath(List(fs))

    def toAbsolutePath(fss: List[FileSystem]): Path = str.toPath(fss).toAbsolutePath
    def toAbsolutePath(fs: FileSystem): Path = str.toAbsolutePath(List(fs))

    def toDirectory(fss: List[FileSystem]): Path = buildPathAsDirectory(str, fss)
    def toDirectory(fs: FileSystem): Path = str.toDirectory(List(fs))

    // TODO this needs to go away because it's gcs specific. Replacing gcs FS with google implementation (when available) will take care of it
    private def buildPathAsDirectory(rawString: String, fileSystems: List[FileSystem]): Path = {
      PathBuilder.findFileSystem(rawString, fileSystems, {
        case fs: GcsFileSystem => Try(fs.getPathAsDirectory(rawString))
        case fs => Try(fs.getPath(rawString))
      })
    }
  }
}
