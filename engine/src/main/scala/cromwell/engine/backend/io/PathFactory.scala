package cromwell.engine.backend.io

import java.nio.file.{Path, FileSystem}

import cromwell.engine.backend.io.filesystem.gcs.GcsFileSystem

import scala.util.{Success, Try}

object PathFactory {
  private def findFileSystem(rawString: String, fss: List[FileSystem], mapping: PartialFunction[FileSystem, Try[Path]]) = {
    fss.toStream collect mapping collectFirst { case Success(p) => p } getOrElse {
      throw new IllegalArgumentException(s"Could not find suitable filesystem to parse $rawString")
    }
  }

  def buildPath(rawString: String, fileSystems: List[FileSystem]): Path = {
    findFileSystem(rawString, fileSystems, {
      case fs: FileSystem => Try(fs.getPath(rawString))
    })
  }

  def buildPathAsDirectory(rawString: String, fileSystems: List[FileSystem]): Path = {
    findFileSystem(rawString, fileSystems, {
      case fs: GcsFileSystem => Try(fs.getPathAsDirectory(rawString))
      case fs => Try(fs.getPath(rawString))
    })
  }
}
