package cloud.nio.spi

import java.nio.ByteBuffer
import java.nio.channels._

class CloudNioWriteChannel(fileProvider: CloudNioFileProvider, retry: CloudNioRetry, cloudNioPath: CloudNioPath)
    extends SeekableByteChannel {
  private var internalPosition: Long = 0
  private val channel: WritableByteChannel = {
    retry.from(
      () => fileProvider.write(cloudNioPath.cloudHost, cloudNioPath.cloudPath)
    )
  }

  override def read(dst: ByteBuffer): Int = throw new NonReadableChannelException

  override def write(src: ByteBuffer): Int = {
    val count = channel.write(src)
    internalPosition += count
    count
  }

  override def position(): Long = internalPosition

  override def position(newPosition: Long): this.type = {
    if (!channel.isOpen)
      throw new ClosedChannelException
    if (newPosition != internalPosition)
      throw new UnsupportedOperationException(s"Cannot change position from $internalPosition to $newPosition")
    this
  }

  override def size(): Long = internalPosition

  override def truncate(size: Long): SeekableByteChannel = throw new UnsupportedOperationException

  override def isOpen: Boolean = channel.isOpen

  override def close(): Unit = channel.close()
}
