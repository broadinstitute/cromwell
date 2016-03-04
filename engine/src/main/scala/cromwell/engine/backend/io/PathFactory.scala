package cromwell.engine.backend.io

import java.nio.file.{Path, FileSystem}

import cromwell.engine.backend.io.filesystem.gcs.GcsFileSystem

import scala.util.{Success, Try}

object PathFactory {

  def buildPath(rawString: String, fileSystems: List[FileSystem]): Path = {
    fileSystems.toStream map { fs => Try(fs.getPath(rawString)) } collectFirst { case Success(p) => p } getOrElse {
      throw new IllegalArgumentException(s"Could not find suitable filesystem to parse $rawString")
    }
  }

  def buildPathAsDirectory(rawString: String, fileSystems: List[FileSystem]): Path = {
    fileSystems.toStream collect {
      case fs: GcsFileSystem => Try(fs.getPathAsDirectory(rawString))
      case fs => Try(fs.getPath(rawString))
    } collectFirst { case Success(p) => p } getOrElse {
      throw new IllegalArgumentException(s"Could not find suitable filesystem to parse $rawString as a Directory")
    }
  }

}
