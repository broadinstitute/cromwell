package cromwell.filesystems.oss.nio

import java.nio.ByteBuffer
import java.nio.channels.{Channels, SeekableByteChannel}

import com.aliyun.oss.OSSClient
import com.aliyun.oss.model.{GenericRequest, GetObjectRequest}

import scala.util.Try

final case class OssFileReadChannel(ossClient: OSSClient, pos: Long, path: OssStoragePath) extends OssFileChannel {
  var internalPosition = pos

  override def position(): Long = {
    synchronized {
      internalPosition
    }
  }

  override def position(newPosition: Long): SeekableByteChannel = {
    if (newPosition < 0) {
      throw new IllegalArgumentException(newPosition.toString)
    }

    synchronized {
      if (newPosition != internalPosition) {
        internalPosition = newPosition
      }

      return this
    }
  }

  override def read(dst: ByteBuffer): Int = {
    dst.mark()

    dst.reset()

    val want = dst.capacity
    val begin: Long = position()
    var end: Long = position + want - 1
    if (begin < 0 || end < 0 || begin > end) {
      throw new IllegalArgumentException(s"being $begin or end $end invalid")
    }

    if (begin >= size) {
        return -1
    }

    if (end >= size()) {
      end = size() - 1
    }

    val getObjectRequest = new GetObjectRequest(path.bucket, path.key)
    getObjectRequest.setRange(begin, end)

    OssStorageRetry.fromTry(
      () => Try{
        val ossObject = ossClient.getObject(getObjectRequest)
        val in = ossObject.getObjectContent
        val channel = Channels.newChannel(in)

        val amt = channel.read(dst)
        channel.close()
        internalPosition += amt
        amt
      }
    )
  }

  override def size(): Long = {
    OssStorageRetry.fromTry(
      () => Try {
        val request = new GenericRequest(path.bucket, path.key)
        request.setLogEnabled(false)
        ossClient.getSimplifiedObjectMeta(request).getSize
      }
    )
  }
}
