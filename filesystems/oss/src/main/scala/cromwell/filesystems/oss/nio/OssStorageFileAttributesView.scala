package cromwell.filesystems.oss.nio

import java.nio.file.NoSuchFileException
import java.nio.file.attribute.{BasicFileAttributeView, FileTime}

import com.aliyun.oss.OSSClient

import scala.util.Try

final case class OssStorageFileAttributesView(ossClient: OSSClient, path: OssStoragePath) extends BasicFileAttributeView {
  override def name(): String = OssStorageFileSystem.URI_SCHEMA

  override def readAttributes(): OssStorageFileAttributes = {
    val ossPath = OssStoragePath.checkPath(path)

    if (ossPath.seemsLikeDirectory) {
      return OssStorageDirectoryAttributes(path)
    }

    if (!ossClient.doesObjectExist(ossPath.bucket, ossPath.key)) {
      throw new NoSuchFileException(path.toString)
    }

    val objectMeta = OssStorageRetry.fromTry(
      () => Try{
        ossClient.getObjectMetadata(path.bucket, path.key)
      }
    )

    OssStorageObjectAttributes(objectMeta, path)
  }

  override def setTimes(lastModifiedTime: FileTime, lastAccessTime: FileTime, createTime: FileTime): Unit = throw new UnsupportedOperationException("OSS object is immutable")
}
