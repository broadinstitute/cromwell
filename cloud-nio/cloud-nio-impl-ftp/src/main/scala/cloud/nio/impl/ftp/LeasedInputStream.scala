package cloud.nio.impl.ftp

import java.io.InputStream

import cats.effect.IO
import cloud.nio.impl.ftp.FtpUtil._
import cloud.nio.impl.ftp.operations.FtpCompletePendingCommand
import io.github.andrebeat.pool.Lease
import org.apache.commons.net.ftp.FTPClient

class LeasedInputStream(cloudHost: String, cloudPath: String, inputStream: InputStream, lease: Lease[FTPClient]) extends InputStream {
  override def read() = inputStream.read()
  override def read(b: Array[Byte]) = inputStream.read(b)
  override def read(b: Array[Byte], off: Int, len: Int): Int = inputStream.read(b, off, len)
  override def skip(n: Long): Long = inputStream.skip(n)
  override def available = inputStream.available()
  override def close() = {
    inputStream.close()
    autoRelease(IO.pure(lease))(FtpCompletePendingCommand(cloudHost, cloudPath, "close input steam").run).void unsafeRunSync()
  }
  override def mark(readlimit: Int) = inputStream.mark(readlimit)
  override def reset() = inputStream.reset()
  override def markSupported = inputStream.markSupported()
}
