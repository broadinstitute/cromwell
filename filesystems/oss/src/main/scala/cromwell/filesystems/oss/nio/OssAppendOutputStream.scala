package cromwell.filesystems.oss.nio

import java.io.{ByteArrayInputStream, OutputStream}

import com.aliyun.oss.OSSClient
import com.aliyun.oss.model.{AppendObjectRequest, GenericRequest}

import scala.util.Try

final case class OssAppendOutputStream(ossClient: OSSClient, path: OssStoragePath, deleteIfExists: Boolean) extends OutputStream {

  var position: Long = {
    val exist = OssStorageRetry.fromTry(
      () => Try{
        val request = new GenericRequest(path.bucket, path.key)
        request.setLogEnabled(false)
        ossClient.doesObjectExist(request)
      }
    )

    var len: Long = 0
    if (exist && deleteIfExists) {
      OssStorageRetry.from(
        () => ossClient.deleteObject(path.bucket, path.key)
      )
    }
    else if (exist) {
      len = OssStorageRetry.from(
        () => ossClient.getObjectMetadata(path.bucket, path.key).getContentLength
      )
    }

    len
  }

  override def write(b: Int): Unit = {
    val arr  = Array[Byte]((b & 0xFF).toByte)

    val appendObjectRequest: AppendObjectRequest = new AppendObjectRequest(path.bucket, path.key, new ByteArrayInputStream(arr))
    this.synchronized {
      appendObjectRequest.setPosition(position)
      val appendObjectResult = OssStorageRetry.fromTry(
        () => Try{
          ossClient.appendObject(appendObjectRequest)
        }
      )

      position = appendObjectResult.getNextPosition()
    }
  }

  override def write(b: Array[Byte]): Unit = {
    val appendObjectRequest: AppendObjectRequest = new AppendObjectRequest(path.bucket, path.key, new ByteArrayInputStream(b))
    this.synchronized {
      appendObjectRequest.setPosition(position)
      val appendObjectResult = OssStorageRetry.fromTry(
        () => Try{
          ossClient.appendObject(appendObjectRequest)
        }
      )
      position = appendObjectResult.getNextPosition()
    }
  }

  override def write(b: Array[Byte], off: Int, len: Int): Unit = {
    if (b == null) {
      throw new NullPointerException
    } else if ((off < 0) || (off > b.length) || (len < 0) || ((off + len) > b.length) || ((off + len) < 0)) {
      throw new IndexOutOfBoundsException
    } else if (len == 0) {
      return
    }

    val s = b.slice(off, off+len)
    val appendObjectRequest: AppendObjectRequest = new AppendObjectRequest(path.bucket, path.key, new ByteArrayInputStream(s))
    this.synchronized {
      appendObjectRequest.setPosition(position)
      val appendObjectResult = OssStorageRetry.fromTry(
        () => Try{
          ossClient.appendObject(appendObjectRequest)
        }
      )
      position = appendObjectResult.getNextPosition()
    }
  }
}
