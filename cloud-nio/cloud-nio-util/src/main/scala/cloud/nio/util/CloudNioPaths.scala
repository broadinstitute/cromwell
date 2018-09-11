package cloud.nio.util

import java.net.{URI, URISyntaxException}
import java.nio.file.{Path, Paths}

import cloud.nio.spi.CloudNioPath

/**
  * Path manipulation methods similar to java.nio.file.Paths.
  *
  * These methods only modify paths but do not connect to the cloud.
  *
  * These methods accept any java.nio.file.Path objects. Unless the Path is a CloudNioPath, the the utilities assume the
  * Path implementations work like a default UnixPath.
  */
object CloudNioPaths {

  /**
    * Parses a path in a way reciprocal with [[cloud.nio.util.CloudNioPaths#showAbsolute(java.nio.file.Path)]].
    *
    * @see [[cloud.nio.util.CloudNioPaths#showAbsolute(java.nio.file.Path)]]
    * @see [[cloud.nio.spi.CloudNioPath#uriAsString()]]
    */
  def get(filePath: String): Path = {
    try {
      // TODO: softer parsing using Guava UrlEscapers. May also be better to list the providers ourselves if possible.
      Paths.get(new URI(filePath))
    } catch {
      case _: URISyntaxException                                               => Paths.get(filePath)
      case iae: IllegalArgumentException if iae.getMessage == "Missing scheme" => Paths.get(filePath)
    }
  }

  /**
    * Return a path in a way reciprocal with [[cloud.nio.util.CloudNioPaths#get]].
    *
    * @see [[cloud.nio.util.CloudNioPaths#get(java.lang.String)]]
    * @see [[cloud.nio.util.CloudNioPaths#showRelative(java.nio.file.Path)]]
    * @see [[cloud.nio.spi.CloudNioPath#uriAsString()]]
    */
  def showAbsolute(path: Path): String = {
    path match {
      case cloudNioPath: CloudNioPath => cloudNioPath.uriAsString
      case _                          => path.toAbsolutePath.toString
    }
  }

  /**
    * When the path is relative returns a relative path in a way reciprocal with resolve.
    * If the path is absolute then it is returned as relative but including the host/bucket.
    *
    * @see [[cloud.nio.util.CloudNioPaths#showAbsolute(java.nio.file.Path)]]
    * @see [[java.nio.file.Path#resolve(java.nio.file.Path)]]
    * @see [[cloud.nio.spi.CloudNioPath#uriAsString()]]
    */
  def showRelative(path: Path): String = {
    path match {
      case cloudNioPath: CloudNioPath => cloudNioPath.relativeDependentPath
      case _ if !path.isAbsolute      => path.normalize().toString
      case _                          => path.getRoot.relativize(path).normalize().toString
    }
  }
}
