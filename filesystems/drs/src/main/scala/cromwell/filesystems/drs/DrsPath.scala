package cromwell.filesystems.drs

import cloud.nio.impl.drs.DrsPathResolver
import cloud.nio.spi.CloudNioPath
import cromwell.core.path.{NioPath, Path}


case class DrsPath(drsPath: CloudNioPath, drsPathResolver: DrsPathResolver) extends Path {

  override protected def nioPath: NioPath = drsPath

  override protected def newPath(nioPath: NioPath): Path = DrsPath(nioPath.asInstanceOf[CloudNioPath], drsPathResolver)

  override def pathAsString: String = drsPath.uriAsString

  override def pathWithoutScheme: String = s"${drsPath.cloudHost}/${drsPath.cloudPath}"
}
