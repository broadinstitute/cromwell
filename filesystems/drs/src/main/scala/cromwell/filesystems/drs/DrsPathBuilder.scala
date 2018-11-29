package cromwell.filesystems.drs

import java.net.URI

import cloud.nio.spi.CloudNioFileSystemProvider
import com.google.common.net.UrlEscapers
import cromwell.core.path.{Path, PathBuilder}

import scala.util.{Failure, Try}


case class DrsPathBuilder(fileSystemProvider: CloudNioFileSystemProvider) extends PathBuilder {

  private val drsScheme: String = fileSystemProvider.getScheme

  override def name: String = "DRS"

  override def build(pathAsString: String): Try[Path] = {
    if (pathAsString.startsWith(s"$drsScheme://")) {
      Try(URI.create(UrlEscapers.urlFragmentEscaper().escape(pathAsString))) flatMap { uri =>
        if (!Option(uri.getScheme).exists(_.equalsIgnoreCase(fileSystemProvider.getScheme))) {
          Failure(new IllegalArgumentException(s"$pathAsString does not have an $drsScheme scheme."))
        } else if (uri.getHost == null && uri.getAuthority == null) {
          Failure(new IllegalArgumentException(s"$pathAsString does not have a valid host."))
        } else {
          Try(DrsPath(fileSystemProvider.getPath(uri)))
        }
      }
    } else {
      Failure(new IllegalArgumentException(s"$pathAsString does not have a $drsScheme} scheme."))
    }
  }
}
