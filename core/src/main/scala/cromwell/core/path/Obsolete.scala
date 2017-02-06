package cromwell.core.path

/**
  * Implements bridges for code that hasn't been updated to use cromwell.core.path.Path.
  *
  * Before adding methods here, recommend checking out and trying [[cromwell.core.path.Path]].
  */
object Obsolete {
  def rm(path: Path) = path.delete(true)

  implicit class PathMethodAliases(val originalPath: Path) extends AnyVal {
    // Instead of `myPath.path`, just use `myPath`.
    def path: Path = originalPath

    // Instead of `myPath.getFileName()`, just use `myPath.getFileName`
    def getFileName(unused: Unit = ()): Path = originalPath.getFileName
  }

  implicit class StringToPath(val path: String) extends AnyVal {
    def toFile: Path = DefaultPathBuilder.get(path)
  }

  object Files {
    def createTempDirectory(prefix: String): DefaultPath = DefaultPathBuilder.createTempDirectory(prefix)

    def exists(path: Path): Boolean = path.exists

    def isSymbolicLink(path: Path): Boolean = path.isSymbolicLink

    // NOTE: When upgrading, symbolicLinkTo returns the target, not the link.
    def createSymbolicLink(link: Path, target: Path): link.type = {
      link.symbolicLinkTo(target)
      link
    }

    def delete(path: Path): Unit = {
      path.delete()
      ()
    }
  }

  object Paths {
    def get(first: String, more: String*): Path = DefaultPathBuilder.get(first, more: _*)
  }

  type File = Path
  val File = ObsoleteFile

  object ObsoleteFile {
    def newTemporaryDirectory(prefix: String = ""): DefaultPath = {
      DefaultPath(better.files.File.newTemporaryDirectory(prefix).path)
    }

    def newTemporaryFile(prefix: String = "", suffix: String = "", parent: Option[Path] = None): Path = {
      parent match {
        case Some(dir) => dir.createTempFile(prefix, suffix)
        case _ => DefaultPathBuilder.createTempFile(prefix, suffix)
      }
    }

    def apply(path: String, fragments: String*) = DefaultPath(better.files.File(path, fragments: _*).path)

    def apply(path: NioPath) = DefaultPath(better.files.File(path).path)

    def apply(path: Path) = path
  }

}
