package cromwell.engine.io

/**
  * Given a list of IoInterfaces, tries to find a suitable one to process the given path.
  */
class IoManager(interfaces: Seq[IoInterface]) extends IoInterface {
  private def findInterface(path: String) = interfaces.find(_.isValidPath(path))
  private def unsupported(path: String) = throw new UnsupportedOperationException(s"No IoInterface able to process $path has been found.")
  private def withInterface[T](path: String, fn: IoInterface => T): T = {
    findInterface(path) map fn getOrElse unsupported(path)
  }

  override def readFile(path: String): String = withInterface(path, _.readFile(path))

  override def writeFile(path: String, content: String): Unit = withInterface(path, _.writeFile(path, content))

  override def listContents(path: String): Iterable[String] = withInterface(path, _.listContents(path))

  override def glob(path: String, pattern: String): Seq[String] = withInterface(path, _.glob(path, pattern))

  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = {
    withInterface(path, _.writeTempFile(path, prefix, suffix, content))
  }

  override def exists(path: String): Boolean = withInterface(path, _.exists(path))

  override def isValidPath(path: String): Boolean = withInterface(path, _.isValidPath(path))

  /*
   * FIXME: use of from is arbitrary here,
   * current implementation does not support copying Local <-> GCS anyway but still not great
   */
  override def copy(from: String, to: String): Unit = withInterface(from, _.copy(from, to))

  override def hash(path: String): String = withInterface(path, _.hash(path))

  override def size(path: String): Long = withInterface(path, _.size(path))
}
