package cloud.nio.impl.ftp

import java.io.OutputStream

import cats.effect.IO
import cloud.nio.impl.ftp.FtpUtil.autoRelease
import cloud.nio.impl.ftp.operations.FtpCompletePendingCommand
import io.github.andrebeat.pool.Lease
import org.apache.commons.net.ftp.FTPClient

class LeasedOutputStream(cloudHost: String, cloudPath: String, outputStream: OutputStream, lease: Lease[FTPClient]) extends OutputStream {
  override def write(b: Int) = outputStream.write(b)
  override def write(b: Array[Byte]) = outputStream.write(b)
  override def write(b: Array[Byte], off: Int, len: Int): Unit = outputStream.write(b, off, len)
  override def flush() = outputStream.flush()
  override def close() = {
    outputStream.close()
    autoRelease(IO.pure(lease))(FtpCompletePendingCommand(cloudHost, cloudPath, "close input steam").run).void unsafeRunSync()
  }
}
