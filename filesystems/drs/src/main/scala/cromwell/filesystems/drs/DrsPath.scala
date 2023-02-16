package cromwell.filesystems.drs

import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import cloud.nio.spi.{CloudNioPath, FileHash}
import cromwell.core.path.{NioPath, Path}

import java.io.IOException


case class DrsPath(drsPath: CloudNioPath, requesterPaysProjectIdOption: Option[String]) extends Path {

  override def nioPath: NioPath = drsPath

  override protected def newPath(nioPath: NioPath): Path = {
    DrsPath(nioPath.asInstanceOf[CloudNioPath], requesterPaysProjectIdOption)
  }

  override def pathAsString: String = drsPath.cloudHost

  override def pathWithoutScheme: String = pathAsString.stripPrefix(drsPath.getFileSystem.provider.getScheme + "://")

  def getFileHash: FileHash = {
    val drsFileSystemProvider = drsPath.getFileSystem.provider.asInstanceOf[DrsCloudNioFileSystemProvider]

    val fileAttributesOption = drsFileSystemProvider.fileProvider.fileAttributes(drsPath.cloudHost, drsPath.cloudPath)

    fileAttributesOption match {
      case Some(fileAttributes) =>
        fileAttributes.fileHash match {
          case Some(fileHash) => fileHash
          case None => throw new IOException(s"Error while resolving DRS path $this. The response from DRS Resolver doesn't contain the 'md5' hash for the file.")
        }
      case None => throw new IOException(s"Error getting file hash of DRS path $this. Reason: File attributes class DrsCloudNioRegularFileAttributes wasn't defined in DrsCloudNioFileProvider.")
    }
  }
}
