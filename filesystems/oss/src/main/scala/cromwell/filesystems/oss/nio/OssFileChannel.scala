package cromwell.filesystems.oss.nio

import java.nio.ByteBuffer
import java.nio.channels.SeekableByteChannel


trait OssFileChannel extends SeekableByteChannel {

  override def isOpen: Boolean = true

  override def close(): Unit = {}

  override def read(dst: ByteBuffer): Int = throw new UnsupportedOperationException()

  override def write(src: ByteBuffer): Int = throw new UnsupportedOperationException()

  override def truncate(size: Long): SeekableByteChannel = throw new UnsupportedOperationException()

}
