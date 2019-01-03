package cloud.nio.impl.ftp

import java.io.{InputStream, OutputStream}

class InputOutputStreams(val inputStream: InputStream, val outputStream: OutputStream) extends AutoCloseable {
  override def close() = {
    try {
      inputStream.close()
    } finally {
      outputStream.close()
    }
  }
}
