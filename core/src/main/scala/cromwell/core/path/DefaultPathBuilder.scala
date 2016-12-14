package cromwell.core.path

import java.net.URI
import java.nio.file.{FileSystems, Path}

import scala.util.Try

/**
  * PathBuilder using the default FileSystem to attempt to build a Path.
  */
case object DefaultPathBuilder extends PathBuilder {
  override def name = "Default"
  override def build(pathAsString: String): Try[Path] = Try {
    val uri = URI.create(pathAsString.replaceAll(" ", "%20"))
    Option(uri.getScheme) match {
      case Some("file") | None => FileSystems.getDefault.getPath(pathAsString)
      case _ => throw new RuntimeException(s"Cannot build a local path from $pathAsString")
    }
  }
}
