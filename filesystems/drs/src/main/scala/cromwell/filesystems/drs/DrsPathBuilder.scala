package cromwell.filesystems.drs

import java.net.URI

import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import com.google.common.net.UrlEscapers
import cromwell.core.path.{Path, PathBuilder}

import scala.util.{Failure, Try}


case class DrsPathBuilder(fileSystemProvider: DrsCloudNioFileSystemProvider) extends PathBuilder {

  private val drsScheme: String = fileSystemProvider.getScheme

  override def name: String = "DRS"

  override def build(pathAsString: String): Try[Path] = {
    if (pathAsString.startsWith(s"$drsScheme://")) {
      Try(URI.create(UrlEscapers.urlFragmentEscaper().escape(pathAsString))) flatMap { uri =>
        if (!Option(uri.getScheme).exists(_.equalsIgnoreCase(fileSystemProvider.getScheme))) {
          Failure(new IllegalArgumentException(s"$pathAsString does not have a $drsScheme scheme."))
        } else if (uri.getHost == null && uri.getAuthority == null) {
          Failure(new IllegalArgumentException(s"$pathAsString does not have a valid host."))
        } else if (uri.getPath == null || uri.getPath.isEmpty || uri.getPath.equalsIgnoreCase("/")) {
          Failure(new IllegalArgumentException(s"$pathAsString does not have a valid path. DRS doesn't support a host only path."))
        } else {
          Try(DrsPath(fileSystemProvider.getPath(uri), fileSystemProvider.drsPathResolver))
        }
      }
    } else {
      Failure(new IllegalArgumentException(s"$pathAsString does not have a $drsScheme scheme."))
    }
  }
}
