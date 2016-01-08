package cromwell.engine.io.gcs

import java.nio.file.Path
import java.nio.file.attribute.{FileTime, BasicFileAttributes}

class GcsFileAttributes(path: Path) extends BasicFileAttributes {
  override def fileKey(): AnyRef = throw new NotImplementedError("To be implemented when/if needed")
  override def isRegularFile: Boolean = throw new NotImplementedError("To be implemented when/if needed")
  override def isOther: Boolean = throw new NotImplementedError("To be implemented when/if needed")
  override def lastModifiedTime(): FileTime = throw new NotImplementedError("To be implemented when/if needed")
  override def size(): Long = throw new NotImplementedError("To be implemented when/if needed")
  override def isDirectory: Boolean = false
  override def isSymbolicLink: Boolean = false
  override def creationTime(): FileTime = throw new NotImplementedError("To be implemented when/if needed")
  override def lastAccessTime(): FileTime = throw new NotImplementedError("To be implemented when/if needed")
}
