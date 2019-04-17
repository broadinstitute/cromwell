package cromwell.filesystems.oss.batch

import cromwell.core.io._
import cromwell.filesystems.oss.OssPath

private case object PartialOssBatchCommandBuilder extends PartialIoCommandBuilder {
  override def sizeCommand = {
    case ossPath: OssPath => OssBatchSizeCommand(ossPath)
  }
  
  override def deleteCommand = {
    case (ossPath: OssPath, swallowIoExceptions) => OssBatchDeleteCommand(ossPath, swallowIoExceptions)
  }
  
  override def copyCommand = {
    case (ossSrc: OssPath, ossDest: OssPath, overwrite) => OssBatchCopyCommand(ossSrc, ossDest, overwrite)
  }
  
  override def hashCommand = {
    case ossPath: OssPath => OssBatchEtagCommand(ossPath)
  }

  override def touchCommand = {
    case ossPath: OssPath => OssBatchTouchCommand(ossPath)
  }

  override def existsCommand = {
    case ossPath: OssPath => OssBatchExistsCommand(ossPath)
  }
}

case object OssBatchCommandBuilder extends IoCommandBuilder(List(PartialOssBatchCommandBuilder))
