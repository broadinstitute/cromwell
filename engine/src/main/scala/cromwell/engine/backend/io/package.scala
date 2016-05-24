package cromwell.engine.backend

import java.nio.file._

import cromwell.engine.backend.io.filesystem.gcs.{ContentTypeOption, GcsFileSystem, GcsFileSystemProvider, NioGcsPath}

import scala.util.{Success, Try}

package object io {
  val defaultGCSFileSystem = GcsFileSystem.defaultGcsFileSystem
  val defaultFileSystem = FileSystems.getDefault
  val defaultFileSystems = List(defaultGCSFileSystem, defaultFileSystem)

  implicit class PathString(val str: String) extends AnyVal {
    def isGcsUrl: Boolean = str.startsWith("gs://")
    def isUriWithProtocol: Boolean = "^[a-z]+://".r.findFirstIn(str).nonEmpty

    def toPath(fss: List[FileSystem]): Path = PathFactory.buildPath(str, fss)
    def toPath(fs: FileSystem): Path = str.toPath(List(fs))

    def toAbsolutePath(fss: List[FileSystem]): Path = str.toPath(fss).toAbsolutePath
    def toAbsolutePath(fs: FileSystem): Path = str.toAbsolutePath(List(fs))

    def toDirectory(fss: List[FileSystem]): Path = PathFactory.buildPathAsDirectory(str, fss)
    def toDirectory(fs: FileSystem): Path = str.toDirectory(List(fs))
  }

  implicit class PathEnhanced(val path: Path) extends AnyVal {
    import better.files._

    def hash = path match {
      case gcs: NioGcsPath => gcs.getFileSystem.provider().asInstanceOf[GcsFileSystemProvider].crc32cHash(gcs)
      case _ => path.md5
    }

    def writeAsJson(content: String): File = {
      Files.write(path, content.getBytes, ContentTypeOption.Json)
    }

    def asDirectory = path.toString.toDirectory(path.getFileSystem)
  }
}
