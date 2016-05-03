package cromwell.engine.backend

import java.nio.file._

import cromwell.core.PathFactory
import cromwell.filesystems.gcs.{NioGcsPath, GcsFileSystemProvider}
import cromwell.filesystems.gcs.{ContentTypeOption, GcsFileSystem}

import scala.util.Try

package object io {
  val defaultGCSFileSystem = GcsFileSystem.defaultGcsFileSystem
  val defaultFileSystem = FileSystems.getDefault
  val defaultFileSystems = List(defaultGCSFileSystem, defaultFileSystem)

  @deprecated("should disappear post PBE")
  private object PathFactoryImpl extends PathFactory
  @deprecated("should disappear post PBE")
  implicit class PathString(val str: String) extends AnyVal {
    def isGcsUrl: Boolean = str.startsWith("gs://")
    def isUriWithProtocol: Boolean = "^[a-z]+://".r.findFirstIn(str).nonEmpty

    def toPath(fss: List[FileSystem]): Path = PathFactoryImpl.buildPath(str, fss)
    def toPath(fs: FileSystem): Path = str.toPath(List(fs))

    def toAbsolutePath(fss: List[FileSystem]): Path = str.toPath(fss).toAbsolutePath
    def toAbsolutePath(fs: FileSystem): Path = str.toAbsolutePath(List(fs))

    def toDirectory(fss: List[FileSystem]): Path = buildPathAsDirectory(str, fss)
    def toDirectory(fs: FileSystem): Path = str.toDirectory(List(fs))

    // TODO this needs to go away because it's gcs specific. Replacing gcs FS with google implementatio (when available) will take care of it
    private def buildPathAsDirectory(rawString: String, fileSystems: List[FileSystem]): Path = {
      PathFactoryImpl.findFileSystem(rawString, fileSystems, {
        case fs: GcsFileSystem => Try(fs.getPathAsDirectory(rawString))
        case fs => Try(fs.getPath(rawString))
      })
    }
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
