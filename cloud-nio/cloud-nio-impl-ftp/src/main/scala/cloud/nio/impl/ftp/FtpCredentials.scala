package cloud.nio.impl.ftp

import java.io.IOException

import org.apache.commons.net.ftp.FTPClient

import scala.util.{Failure, Success, Try}

sealed trait FtpCredentials {
  def login(ftpClient: FTPClient) = {}
}

// Yes, FTP uses plain text username / password
case class FtpAuthenticatedCredentials(username: String, password: String, account: Option[String]) extends FtpCredentials {
  override def login(ftpClient: FTPClient) = {
    lazy val replyString = Option(ftpClient.getReplyString).getOrElse("N/A")

    def loginAction = account match {
      case Some(a) => ftpClient.login(username, password, a)
      case _ => ftpClient.login(username, password)
    }

    Try(loginAction) match {
      case Success(true) =>
      case Success(false) =>
        throw new IOException(s"Failed to authenticate to FTP server: $replyString")
      case Failure(f) =>
        throw new IOException(s"Failed to authenticate to FTP server: $replyString", f)
    }
  }
}

case object FtpAnonymousCredentials extends FtpCredentials
