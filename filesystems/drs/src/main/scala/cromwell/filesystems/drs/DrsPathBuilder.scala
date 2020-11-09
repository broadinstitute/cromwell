package cromwell.filesystems.drs

import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import cromwell.core.path.{Path, PathBuilder}

import scala.util.{Failure, Try}


case class DrsPathBuilder(fileSystemProvider: DrsCloudNioFileSystemProvider,
                          requesterPaysProjectIdOption: Option[String]) extends PathBuilder {

  private val drsScheme: String = fileSystemProvider.getScheme

  override def name: String = "DRS"

  override def build(pathAsString: String): Try[Path] = {
    if (pathAsString.startsWith(s"$drsScheme://")) {
      Try(DrsPath(fileSystemProvider.getCloudNioPath(pathAsString), requesterPaysProjectIdOption))
    } else {
      Failure(new IllegalArgumentException(s"$pathAsString does not have a $drsScheme scheme."))
    }
  }
}
