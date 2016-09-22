package cromwell.core.path

import java.nio.file.Path

object PathImplicits {
  implicit class EnhancedPath(val path: Path) extends AnyVal {
    def swapExt(oldExt: String, newExt: String): Path = {
      path.getFileSystem.getPath(s"${path.toString.stripSuffix(oldExt)}$newExt")
    }

    def untailed = UntailedWriter(path)

    def tailed(tailedSize: Int) = TailedWriter(path, tailedSize)
  }
}
