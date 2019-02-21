package cromwell.filesystems.sra

import cromwell.core.path.{NioPath, Path, PathBuilder}

import scala.util.{Failure, Success, Try}

class SraPathBuilder extends PathBuilder {
  override def name: String = "NIH Sequence Read Archive"

  def build(path: String): Try[SraPath] = {
    if (!path.startsWith(SraPath.Scheme)) {
      return Failure(new IllegalArgumentException(s"$path is not a SRA path"))
    }
    val parts = path.stripPrefix(SraPath.Scheme).split("/", 2)
    if (parts.length != 2) {
      return Failure(new IllegalArgumentException(s"$path does not have enough components"))
    }
    Success(new SraPath(parts(0), parts(1)))
  }
}

object SraPath {
  val Scheme = "sra://"
}

case class SraPath(accession: String, path: String) extends Path {
  override def pathAsString: String = SraPath.Scheme + pathWithoutScheme
  override def pathWithoutScheme: String = accession + "/" + path

  protected def newPath(nioPath: NioPath): Path = throw new UnsupportedOperationException("'newPath' not implemented for SraPath")
  protected def nioPath: NioPath = throw new UnsupportedOperationException("'nioPath' not implemented for SraPath")
}
