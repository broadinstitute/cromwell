package cloud.nio.util

import java.io.{ByteArrayOutputStream, InputStream}
import java.nio.ByteBuffer
import java.nio.channels.{Channels, Pipe, WritableByteChannel}

import scala.util.Try

/**
  * Various wrappers for creating ByteChannels on background threads.
  *
  * Sounds a lot like Netty.
  * - https://www.safaribooksonline.com/library/view/netty-in-action/9781617291470/
  * - https://netty.io/4.1/api/io/netty/buffer/Unpooled.html
  * - http://www.baeldung.com/netty
  */
object ChannelUtil {

  /**
    * Creates a WritableByteChannel from an InputStream updated on a background thread.
    *
    * @param threadName The name of the background thread
    * @param consumer   Consumes the input stream.
    * @return A writable byte channel for the buffer
    */
  def pipedStreamWriter(threadName: String)(consumer: InputStream => Unit): WritableByteChannel = {
    val pipe = Pipe.open()
    var threadResult: Option[Try[Unit]] = None
    val runnable: Runnable = () => {
      threadResult = Option(
        Try(
          consumer(Channels.newInputStream(pipe.source))
        )
      )
    }
    val thread = new Thread(runnable, threadName)
    thread.setDaemon(true)
    thread.start()

    new WritableByteChannel {
      override def write(src: ByteBuffer): Int = pipe.sink.write(src)

      override def isOpen: Boolean = pipe.sink.isOpen

      override def close(): Unit = {
        pipe.sink.close()
        thread.join()
        pipe.source.close()
        val _ = threadResult.get.get
      }
    }

  }

  /**
    * Creates an unlimited buffer of memory to consume data before passing it into the consumer.
    *
    * Useful only for uploaders that must know the full size ahead of time.
    *
    * @param threadName The name of the background thread
    * @param consumer   Consumes the buffered data
    * @return A writable byte channel for the buffer
    */
  def blockingMemoryWriter(threadName: String)(consumer: Array[Byte] => Unit): WritableByteChannel = {
    var threadResult: Option[Try[Unit]] = None
    val byteArrayOutputStream = new ByteArrayOutputStream
    var done = false
    val bufferLock = new Object
    val runnable: Runnable = () => {
      bufferLock.synchronized {
        if (!done)
          bufferLock.wait()
      }
      threadResult = Option(
        Try(
          consumer(byteArrayOutputStream.toByteArray)
        )
      )
    }

    val thread = new Thread(runnable, threadName)
    thread.start()

    val wrapper = Channels.newChannel(byteArrayOutputStream)

    new WritableByteChannel {
      override def write(src: ByteBuffer): Int = wrapper.write(src)

      override def isOpen: Boolean = wrapper.isOpen

      override def close(): Unit = {
        wrapper.close()
        bufferLock.synchronized {
          done = true
          bufferLock.notify()
        }
        thread.join()
        val _ = threadResult.get.get
      }
    }

  }

}
