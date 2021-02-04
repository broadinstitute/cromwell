package cromwell.filesystems.drs

import cloud.nio.spi.CloudNioPath
import cromwell.core.path.{NioPath, Path}


case class DrsPath(drsPath: CloudNioPath, requesterPaysProjectIdOption: Option[String]) extends Path {

  override protected def nioPath: NioPath = drsPath

  override protected def newPath(nioPath: NioPath): Path = {
    DrsPath(nioPath.asInstanceOf[CloudNioPath], requesterPaysProjectIdOption)
  }

  override def pathAsString: String = drsPath.cloudHost

  override def pathWithoutScheme: String = pathAsString.stripPrefix(drsPath.getFileSystem.provider.getScheme + "://")
}
