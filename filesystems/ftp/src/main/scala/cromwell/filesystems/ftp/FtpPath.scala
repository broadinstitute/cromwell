package cromwell.filesystems.ftp

import cloud.nio.spi.CloudNioPath
import cromwell.core.path.{NioPath, Path}

case class FtpPath(ftpPath: CloudNioPath) extends Path {
  override protected def nioPath = ftpPath
  override protected def newPath(nioPath: NioPath) = FtpPath(nioPath.asInstanceOf[CloudNioPath])
  override def pathAsString = nioPath.uriAsString
  override def pathWithoutScheme = nioPath.cloudPath
  override def createPermissionedDirectories() = super.createDirectories()
}
