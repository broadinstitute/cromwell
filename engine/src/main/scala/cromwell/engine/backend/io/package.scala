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
}
