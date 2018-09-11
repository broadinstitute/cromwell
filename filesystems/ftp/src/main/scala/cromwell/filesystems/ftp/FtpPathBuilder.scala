package cromwell.filesystems.ftp

import java.net.URI

import cloud.nio.spi.CloudNioFileSystemProvider
import com.google.common.net.UrlEscapers
import cromwell.core.path.PathBuilder
import org.slf4j.LoggerFactory

import scala.util.{Failure, Try}

object FtpPathBuilder {
  val logger = LoggerFactory.getLogger("FTPLogger")
}

case class FtpPathBuilder(fileSystemProvider: CloudNioFileSystemProvider) extends PathBuilder {
  override def name = "FTP"

  override def build(string: String) = {
    if (string == "ftp://") Failure(new IllegalArgumentException(s"$string does not have a valid host"))
    else {
      Try(URI.create(UrlEscapers.urlFragmentEscaper().escape(string))) flatMap { uri =>
        if (!Option(uri.getScheme).exists(_.equalsIgnoreCase(fileSystemProvider.getScheme))) {
          Failure(new IllegalArgumentException(s"$string does not have an ftp scheme"))
        } else if (uri.getHost == null) {
          Failure(new IllegalArgumentException(s"$string does not have a valid host"))
        } else {
          Try {
            FtpPath(fileSystemProvider.getPath(uri))
          }
        }
      }
    }
  }
}
