package cromwell.filesystems.oss.nio

import java.nio.file.NoSuchFileException
import java.nio.file.attribute.{BasicFileAttributeView, FileTime}

import com.aliyun.oss.OSSClient
import com.aliyun.oss.model.GenericRequest

import scala.util.Try

final case class OssStorageFileAttributesView(ossClient: OSSClient, path: OssStoragePath) extends BasicFileAttributeView {
  override def name(): String = OssStorageFileSystem.URI_SCHEMA

  override def readAttributes(): OssStorageFileAttributes = {
    val ossPath = OssStoragePath.checkPath(path)

    if (ossPath.seemsLikeDirectory) {
      return OssStorageDirectoryAttributes(path)
    }

    val request = new GenericRequest(ossPath.bucket, ossPath.key)
    request.setLogEnabled(false)
    if (!ossClient.doesObjectExist(request)) {
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
