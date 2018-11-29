package cromwell.filesystems.drs

import cloud.nio.spi.CloudNioPath
import cromwell.core.path.{NioPath, Path}


case class DrsPath(drsPath: CloudNioPath) extends Path {

  override protected def nioPath: NioPath = drsPath

  override protected def newPath(nioPath: NioPath): Path = DrsPath(nioPath.asInstanceOf[CloudNioPath])

  override def pathAsString: String = drsPath.uriAsString

  override def pathWithoutScheme: String = drsPath.cloudPath
}
