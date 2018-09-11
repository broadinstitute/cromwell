package cloud.nio.util

import java.nio.file.{Files, Path}

import scala.collection.JavaConverters._

/**
  * Common file access utilities similar to java.nio.file.Files.
  */
object CloudNioFiles {

  /**
    * Lists all files under a path.
    */
  def listRegularFiles(path: Path): Iterator[Path] = {
    Files
      .walk(path, Int.MaxValue)
      .iterator
      .asScala
      .filter(Files.isRegularFile(_))
  }

  /**
    * Returns an iterator of all regular files under sourcePath mapped relatively to targetPath.
    */
  def relativeFiles(
    sourcePath: Path,
    targetPath: Path,
    pathToString: Path => String = CloudNioPaths.showRelative
  ): Iterator[(Path, Path)] = {
    def relativeSourceTarget(path: Path): (Path, Path) = {
      val sourceRelative = sourcePath.relativize(path)
      val relativeString = pathToString(sourceRelative)
      val targetResolve = targetPath.resolve(relativeString)
      (path, targetResolve)
    }

    listRegularFiles(sourcePath).map(relativeSourceTarget)
  }

}
