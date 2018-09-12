package cloud.nio.spi

import java.io.FileNotFoundException
import java.nio.ByteBuffer
import java.nio.channels.{
  ClosedChannelException,
  NonWritableChannelException,
  ReadableByteChannel,
  SeekableByteChannel
}

class CloudNioReadChannel(fileProvider: CloudNioFileProvider, retry: CloudNioRetry, cloudNioPath: CloudNioPath)
  extends SeekableByteChannel {
  private var internalPosition: Long = 0
  private var channel: ReadableByteChannel = channelPosition(0)

  override def read(dst: ByteBuffer): Int = {
    var resetConnection = false
    val count = retry.from(
      () => {
        try {
          if (resetConnection) {
            if (channel.isOpen) channel.close()
            channel = fileProvider.read(cloudNioPath.cloudHost, cloudNioPath.cloudPath, internalPosition)
          }
          channel.read(dst)
        } catch {
          case exception: Exception =>
            resetConnection = true
            throw exception
        }
      }
    )
    if (count > 0)
      internalPosition += count
    count
  }

  override def write(src: ByteBuffer): Int = throw new NonWritableChannelException

  override def position(): Long = internalPosition

  override def position(newPosition: Long): this.type = {
    if (!channel.isOpen)
      throw new ClosedChannelException

    channel.close()
    internalPosition = newPosition
    channel = channelPosition(newPosition)
    this
  }

  private def channelPosition(newPosition: Long): ReadableByteChannel = {
    retry.from(
      () => fileProvider.read(cloudNioPath.cloudHost, cloudNioPath.cloudPath, newPosition)
    )
  }

  override def size(): Long = {
    retry
      .from(
        () => fileSize
      )
      .getOrElse(throw new FileNotFoundException(cloudNioPath.uriAsString))
  }

  override def truncate(size: Long): SeekableByteChannel = throw new NonWritableChannelException

  override def isOpen: Boolean = channel.isOpen

  override def close(): Unit = channel.close()
  
  protected def fileSize = fileProvider.fileAttributes(cloudNioPath.cloudHost, cloudNioPath.cloudPath).map(_.size())
}
