package cromwell.core.path

import java.nio.file.Path

object PathImplicits {
  implicit class EnhancedPath(val path: Path) extends AnyVal {
    def swapExt(oldExt: String, newExt: String): Path = {
      path.getFileSystem.getPath(s"${path.toString.stripSuffix(oldExt)}$newExt")
    }

    def untailed = UntailedWriter(path)

    def tailed(tailedSize: Int) = TailedWriter(path, tailedSize)

    def toRealString: String = java.net.URLDecoder.decode(path.toUri.toString, "UTF-8")
  }
}
