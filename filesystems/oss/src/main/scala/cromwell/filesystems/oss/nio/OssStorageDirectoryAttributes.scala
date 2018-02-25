package cromwell.filesystems.oss.nio

import java.nio.file.attribute.FileTime

import scala.collection.mutable.Map

final case class OssStorageDirectoryAttributes(path: OssStoragePath) extends OssStorageFileAttributes {
  override def creationTime(): FileTime = FileTime.fromMillis(0)

  override def lastAccessTime(): FileTime = FileTime.fromMillis(0)

  override def lastModifiedTime(): FileTime = creationTime()

  override def isRegularFile: Boolean = false

  override def isDirectory: Boolean = true

  override def isSymbolicLink: Boolean = false

  override def isOther: Boolean = false

  override def size(): Long = 0

  override def fileKey(): AnyRef = path.pathAsString

  override def expires: FileTime = FileTime.fromMillis(0)

  override def cacheControl(): Option[String] = None

  override def contentDisposition: Option[String] = None

  override def contentEncoding: Option[String] = None

  override def etag: Option[String] = None

  override def userMeta: Map[String, String] = Map.empty[String, String]

}
