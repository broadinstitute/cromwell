package cromwell.core.path

import java.nio.file.FileSystems

import scala.util.Try

class MappedPathBuilder(prefix: String, mappedRoot: String) extends PathBuilder {
  override def name: String = s"$prefix mapped to $mappedRoot"

  override def build(pathAsString: String): Try[Path] = Try {
    if (pathAsString.startsWith(prefix)) {
      val fileSystem = FileSystems.getDefault
      val root = fileSystem.getPath(mappedRoot)
      val path = fileSystem.getPath(mappedRoot, pathAsString.stripPrefix(prefix))
      MappedPath(prefix, root, path)

    } else {
      throw new IllegalArgumentException(s"$pathAsString must start with $prefix")
    }
  }
}

case class MappedPath(prefix: String, mappedRoot: NioPath, nioPath: NioPath) extends Path {
  override protected def newPath(nioPath: NioPath): Path = MappedPath(prefix, mappedRoot, nioPath)

  override def pathAsString: String = nioPath.toString

  override def pathWithoutScheme: String = mappedRoot.relativize(nioPath).toString

  def prefixedPathAsString: String = prefix + pathWithoutScheme
}
