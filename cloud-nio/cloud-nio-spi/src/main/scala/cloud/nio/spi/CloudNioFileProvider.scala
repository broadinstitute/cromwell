package cloud.nio.spi

import java.nio.channels.{ReadableByteChannel, WritableByteChannel}
import java.nio.file.attribute.FileTime

abstract class CloudNioFileProvider {
  def existsPath(cloudHost: String, cloudPath: String): Boolean

  def existsPaths(cloudHost: String, cloudPathPrefix: String): Boolean

  /**
    * Returns a listing of keys within a bucket starting with prefix. The returned keys should include the prefix. The
    * paths must be absolute, but the key should not begin with a slash.
    */
  def listObjects(cloudHost: String, cloudPathPrefix: String, markerOption: Option[String]): CloudNioFileList

  def copy(sourceCloudHost: String, sourceCloudPath: String, targetCloudHost: String, targetCloudPath: String): Unit

  def deleteIfExists(cloudHost: String, cloudPath: String): Boolean

  def read(cloudHost: String, cloudPath: String, offset: Long): ReadableByteChannel

  def write(cloudHost: String, cloudPath: String): WritableByteChannel

  def fileAttributes(cloudHost: String, cloudPath: String): Option[CloudNioRegularFileAttributes]

  def createDirectory(cloudHost: String, cloudPath: String): Unit = {}
}

object CloudNioFileProvider {
  val UnknownTime: FileTime = FileTime.fromMillis(0L)
  val UnknownSize: Long = 0L
}
