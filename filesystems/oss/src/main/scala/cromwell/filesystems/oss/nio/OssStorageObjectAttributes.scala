package cromwell.filesystems.oss.nio

import java.nio.file.attribute.FileTime

import com.aliyun.oss.model.ObjectMetadata

import scala.jdk.CollectionConverters._
import scala.collection.mutable.Map
import scala.util.Try

final case class OssStorageObjectAttributes(objectMeta: ObjectMetadata, path: OssStoragePath) extends OssStorageFileAttributes {
  override def creationTime(): FileTime = {
    FileTime.fromMillis(objectMeta.getLastModified.getTime)
  }

  override def lastAccessTime(): FileTime = FileTime.fromMillis(0)

  override def lastModifiedTime(): FileTime = creationTime()

  override def isRegularFile: Boolean = true

  override def isDirectory: Boolean = false

  override def isSymbolicLink: Boolean = false

  override def isOther: Boolean = false

  override def size(): Long = objectMeta.getContentLength

  override def fileKey(): AnyRef = path.pathAsString

  // oss sdk has an issule: throw NullPointerException when no expire time exists.
  override def expires: FileTime = FileTime.fromMillis(Try{objectMeta.getExpirationTime.getTime} getOrElse (0))

  override def cacheControl(): Option[String] = Option(objectMeta.getCacheControl)

  override def contentDisposition: Option[String] = Option(objectMeta.getContentDisposition)

  override def contentEncoding: Option[String] = Option(objectMeta.getContentEncoding)

  override def etag: Option[String] = Option(objectMeta.getETag)

  override def userMeta: Map[String, String] = objectMeta.getUserMetadata.asScala
}
