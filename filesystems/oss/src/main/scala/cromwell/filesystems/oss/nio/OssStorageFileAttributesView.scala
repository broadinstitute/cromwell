package cromwell.filesystems.oss.nio

import java.nio.file.NoSuchFileException
import java.nio.file.attribute.{BasicFileAttributeView, FileTime}
import java.util.Date

import com.aliyun.oss.OSSClient
import com.aliyun.oss.model.{CopyObjectRequest, CopyObjectResult, GenericRequest}

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

  override def setTimes(lastModifiedTime: FileTime, lastAccessTime: FileTime, createTime: FileTime): Unit = {
    val meta = ossClient.getObjectMetadata(path.bucket, path.key)
    meta.setLastModified(new Date(lastModifiedTime.toMillis))

    val copyReq = new CopyObjectRequest(path.bucket, path.key, path.bucket, path.key)
    copyReq.setNewObjectMetadata(meta)
    val _: CopyObjectResult = ossClient.copyObject(copyReq)
  }
}
