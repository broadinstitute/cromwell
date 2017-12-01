package cromwell.filesystems.gcs.batch

import cromwell.core.io._
import cromwell.core.path.BetterFileMethods.OpenOptions
import cromwell.core.path.Path
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.oss.OssPath
import cromwell.filesystems.oss._

trait GcsBatchCommandBuilder extends DefaultIoCommandBuilder {
  override def sizeCommand(path: Path) = path match {
    case ossPath: OssPath => OssIoSizeCommand(ossPath)
    case gcsPath: GcsPath => GcsBatchSizeCommand(gcsPath)
    case _ => super.sizeCommand(path)
  }
  
  override def deleteCommand(path: Path, swallowIoExceptions: Boolean = false) = path match {
    case ossPath: OssPath => OssIoDeleteCommand(ossPath, swallowIoExceptions)
    case gcsPath: GcsPath => GcsBatchDeleteCommand(gcsPath, swallowIoExceptions)
    case _ => super.deleteCommand(path, swallowIoExceptions)
  }
  
  override def copyCommand(src: Path, dest: Path, overwrite: Boolean = true) =  (src, dest) match {
    case (ossSrc: OssPath, ossDest: OssPath) => OssIoCopyCommand(ossSrc, ossDest, overwrite)
    case (gcsSrc: GcsPath, gcsDest: GcsPath) => GcsBatchCopyCommand(gcsSrc, gcsDest, overwrite)
    case _ => super.copyCommand(src, dest, overwrite)
  }
  
  override def hashCommand(path: Path) =  path match {
    case ossPath: OssPath => OssIoHashCommand(ossPath)
    case gcsPath: GcsPath => GcsBatchCrc32Command(gcsPath)
    case _ => super.hashCommand(path)
  }

  override def touchCommand(path: Path) =  path match {
    case ossPath: OssPath => OssIoTouchCommand(ossPath)
    case gcsPath: GcsPath => GcsBatchTouchCommand(gcsPath)
    case _ => super.touchCommand(path)
  }

  override def writeCommand(path: Path, content: String, options: OpenOptions) = path match {
    case ossPath: OssPath => OssIoWriteCommand(ossPath, content, options)
    case _ => super.writeCommand(path, content, options)
  }

  override def contentAsStringCommand(path: Path) = path match {
    case ossPath: OssPath => OssIoContentAsStringCommand(ossPath)
    case _ => super.contentAsStringCommand(path)
  }
}
