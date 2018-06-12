package cromwell.filesystems.gcs.batch

import cromwell.core.io._
import cromwell.filesystems.gcs.GcsPath

private case object PartialGcsBatchCommandBuilder extends PartialIoCommandBuilder {
  override def sizeCommand = {
    case gcsPath: GcsPath => GcsBatchSizeCommand(gcsPath)
  }
  
  override def deleteCommand = {
    case (gcsPath: GcsPath, swallowIoExceptions) => GcsBatchDeleteCommand(gcsPath, swallowIoExceptions)
  }
  
  override def copyCommand = {
    case (gcsSrc: GcsPath, gcsDest: GcsPath, overwrite) => GcsBatchCopyCommand(gcsSrc, gcsDest, overwrite)
  }
  
  override def hashCommand = {
    case gcsPath: GcsPath => GcsBatchCrc32Command(gcsPath)
  }

  override def touchCommand = {
    case gcsPath: GcsPath => GcsBatchTouchCommand(gcsPath)
  }

  override def existsCommand = {
    case gcsPath: GcsPath => GcsBatchExistsCommand(gcsPath)
  }

  override def isDirectoryCommand = {
    case gcsPath: GcsPath => GcsBatchIsDirectoryCommand(gcsPath)
  }
}

case object GcsBatchCommandBuilder extends IoCommandBuilder(List(PartialGcsBatchCommandBuilder))
