package cromwell.engine.io.gcs

import java.nio.file.Path
import java.nio.file.attribute.{FileTime, BasicFileAttributes}

class GcsFileAttributes(path: Path) extends BasicFileAttributes {

  override def fileKey(): AnyRef = ???

  override def isRegularFile: Boolean = ???

  override def isOther: Boolean = ???

  override def lastModifiedTime(): FileTime = ???

  override def size(): Long = ???

  override def isDirectory: Boolean = false

  override def isSymbolicLink: Boolean = false

  override def creationTime(): FileTime = ???

  override def lastAccessTime(): FileTime = ???
}
