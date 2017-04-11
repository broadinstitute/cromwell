package cromwell.filesystems.gcs.batch

import cromwell.core.io._
import cromwell.core.path.Path
import cromwell.filesystems.gcs.GcsPath

trait GcsBatchCommandBuilder extends DefaultIoCommandBuilder {
  override def sizeCommand(path: Path) = path match {
    case gcsPath: GcsPath => GcsBatchSizeCommand(gcsPath)
    case _ => super.sizeCommand(path)
  }
  
  override def deleteCommand(path: Path, swallowIoExceptions: Boolean = false) = path match {
    case gcsPath: GcsPath => GcsBatchDeleteCommand(gcsPath, swallowIoExceptions)
    case _ => super.deleteCommand(path, swallowIoExceptions)
  }
  
  override def copyCommand(src: Path, dest: Path, overwrite: Boolean = true) =  (src, dest) match {
    case (gcsSrc: GcsPath, gcsDest: GcsPath) => GcsBatchCopyCommand(gcsSrc, gcsDest, overwrite)
    case _ => super.copyCommand(src, dest, overwrite)
  }
  
  override def hashCommand(path: Path) =  path match {
    case gcsPath: GcsPath => GcsBatchCrc32Command(gcsPath)
    case _ => super.hashCommand(path)
  }
}
