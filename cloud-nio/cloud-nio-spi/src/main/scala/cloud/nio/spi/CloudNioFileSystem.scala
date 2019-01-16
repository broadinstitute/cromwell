package cloud.nio.spi

import java.nio.file._
import java.nio.file.attribute.UserPrincipalLookupService

import scala.collection.JavaConverters._

object CloudNioFileSystem {
  val Separator: String = "/"
}

//TODO: Use stronger type for CloudNioFileSystemProvider to avoid casting
class CloudNioFileSystem(override val provider: CloudNioFileSystemProvider, val host: String) extends FileSystem {

  override def getPath(first: String, more: String*): CloudNioPath = getCloudNioPath(UnixPath.getPath(first, more: _*))

  private[this] def getCloudNioPath(unixPath: UnixPath): CloudNioPath = new CloudNioPath(this, unixPath)

  override def close(): Unit = {
    // do nothing currently.
  }

  override def isOpen: Boolean = {
    true
  }

  override def isReadOnly: Boolean = {
    false
  }

  override def getSeparator: String = {
    CloudNioFileSystem.Separator
  }

  override def getRootDirectories: java.lang.Iterable[Path] = {
    Set[Path](getPath(UnixPath.Root)).asJava
  }

  override def getFileStores: java.lang.Iterable[FileStore] = {
    Set.empty[FileStore].asJava
  }

  override def getPathMatcher(syntaxAndPattern: String): PathMatcher = {
    FileSystems.getDefault.getPathMatcher(syntaxAndPattern)
  }

  override def getUserPrincipalLookupService: UserPrincipalLookupService = {
    throw new UnsupportedOperationException
  }

  override def newWatchService(): WatchService = {
    throw new UnsupportedOperationException
  }

  override def supportedFileAttributeViews(): java.util.Set[String] = {
    Set("basic", CloudNioFileAttributeView.Name).asJava
  }

  def canEqual(other: Any): Boolean = other.isInstanceOf[CloudNioFileSystem]

  override def equals(other: Any): Boolean = other match {
    case that: CloudNioFileSystem =>
      (that canEqual this) &&
        provider == that.provider &&
        host == that.host
    case _ => false
  }

  override def hashCode(): Int = {
    val state = List(provider, host)
    state.map(_.hashCode()).foldLeft(0)((a, b) => 31 * a + b)
  }
}
