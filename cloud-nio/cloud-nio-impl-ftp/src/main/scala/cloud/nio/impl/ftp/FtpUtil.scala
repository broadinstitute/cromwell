package cloud.nio.impl.ftp

import java.io.IOException

import cats.effect.{ExitCase, IO}
import io.github.andrebeat.pool.Lease
import org.apache.commons.net.ftp.{FTPClient, FTPReply}

object FtpUtil {
  implicit class EnhancedCloudPath(val cloudPath: String) extends AnyVal {
    def ensureSlashedPrefix = if (cloudPath.startsWith("/")) cloudPath else s"/$cloudPath"
    def ensureSlashed = if (cloudPath.endsWith("/")) cloudPath else s"$cloudPath/"
  }

  case class FtpIoException(message: String, code: Int, replyString: String, cause: Option[Throwable] = None)
      extends IOException(s"$message: $replyString", cause.orNull) {
    def isTransient = FTPReply.isNegativeTransient(code)
    def isFatal = FTPReply.isNegativePermanent(code)
  }

  def autoRelease[A](acquire: IO[Lease[FTPClient]])(action: FTPClient => IO[A]): IO[A] =
    acquire.bracketCase(lease => action(lease.get())) {
      // If there's a cause, the call to the FTP client threw an exception, assume the connection is compromised and invalidate the lease
      case (lease, ExitCase.Error(FtpIoException(_, _, _, Some(_)))) => IO(lease.invalidate())
      case (lease, _) => IO(lease.release())
    }
}
