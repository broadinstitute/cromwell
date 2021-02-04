package cromwell.filesystems.gcs.batch

import cromwell.core.io._
import cromwell.core.path.Path
import cromwell.filesystems.gcs.GcsPath

import scala.util.Try

private case object PartialGcsBatchCommandBuilder extends PartialIoCommandBuilder {
  override def sizeCommand: PartialFunction[Path, Try[GcsBatchSizeCommand]] = {
    case gcsPath: GcsPath => GcsBatchSizeCommand.forPath(gcsPath)
  }
  
  override def deleteCommand: PartialFunction[(Path, Boolean), Try[GcsBatchDeleteCommand]] = {
    case (gcsPath: GcsPath, swallowIoExceptions) => GcsBatchDeleteCommand.forPath(gcsPath, swallowIoExceptions)
  }
  
  override def copyCommand: PartialFunction[(Path, Path), Try[GcsBatchCopyCommand]] = {
    case (gcsSrc: GcsPath, gcsDest: GcsPath) => GcsBatchCopyCommand.forPaths(gcsSrc, gcsDest)
  }
  
  override def hashCommand: PartialFunction[Path, Try[GcsBatchCrc32Command]] = {
    case gcsPath: GcsPath => GcsBatchCrc32Command.forPath(gcsPath)
  }

  override def touchCommand: PartialFunction[Path, Try[GcsBatchTouchCommand]] = {
    case gcsPath: GcsPath => GcsBatchTouchCommand.forPath(gcsPath)
  }

  override def existsCommand: PartialFunction[Path, Try[GcsBatchExistsCommand]] = {
    case gcsPath: GcsPath => GcsBatchExistsCommand.forPath(gcsPath)
  }

  override def isDirectoryCommand: PartialFunction[Path, Try[GcsBatchIsDirectoryCommand]] = {
    case gcsPath: GcsPath => GcsBatchIsDirectoryCommand.forPath(gcsPath)
  }
}

case object GcsBatchCommandBuilder extends IoCommandBuilder(List(PartialGcsBatchCommandBuilder))
