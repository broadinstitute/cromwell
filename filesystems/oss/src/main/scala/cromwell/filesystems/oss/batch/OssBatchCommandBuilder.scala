package cromwell.filesystems.oss.batch

import cromwell.core.io._
import cromwell.core.path.Path
import cromwell.filesystems.oss.OssPath

import scala.util.Try

private case object PartialOssBatchCommandBuilder extends PartialIoCommandBuilder {
  override def sizeCommand: PartialFunction[Path, Try[IoSizeCommand]] = {
    case ossPath: OssPath => Try(OssBatchSizeCommand(ossPath))
  }
  
  override def deleteCommand: PartialFunction[(Path, Boolean), Try[IoDeleteCommand]] = {
    case (ossPath: OssPath, swallowIoExceptions) => Try(OssBatchDeleteCommand(ossPath, swallowIoExceptions))
  }
  
  override def copyCommand: PartialFunction[(Path, Path), Try[IoCopyCommand]] = {
    case (ossSrc: OssPath, ossDest: OssPath) => Try(OssBatchCopyCommand(ossSrc, ossDest))
  }
  
  override def hashCommand: PartialFunction[Path, Try[IoHashCommand]] = {
    case ossPath: OssPath => Try(OssBatchEtagCommand(ossPath))
  }

  override def touchCommand: PartialFunction[Path, Try[IoTouchCommand]] = {
    case ossPath: OssPath => Try(OssBatchTouchCommand(ossPath))
  }

  override def existsCommand: PartialFunction[Path, Try[IoExistsCommand]] = {
    case ossPath: OssPath => Try(OssBatchExistsCommand(ossPath))
  }
}

case object OssBatchCommandBuilder extends IoCommandBuilder(List(PartialOssBatchCommandBuilder))
