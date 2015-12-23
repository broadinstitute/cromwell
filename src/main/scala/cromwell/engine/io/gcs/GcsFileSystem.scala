package cromwell.engine.io.gcs

import java.lang.Iterable
import java.nio.file._
import java.nio.file.attribute.UserPrincipalLookupService
import java.nio.file.spi.FileSystemProvider
import java.util

import scala.collection.JavaConverters._
import scala.collection.immutable.HashSet
import scala.concurrent.ExecutionContext
import scala.util.Try

object GcsFileSystem {
  def getInstance(provider: GcsFileSystemProvider) = new GcsFileSystem(provider)
  def getInstance(interface: GoogleCloudStorage)(implicit executionContext: ExecutionContext) = new GcsFileSystem(GcsFileSystemProvider.getInstance(interface))
}

class GcsFileSystem(gcsFileSystemProvider: GcsFileSystemProvider) extends FileSystem {

  private val SEPARATOR = "/"
  private [io] val PROTOCOL = "gs://"
  private val gsUriRegex = s"""$PROTOCOL(.*)""".r

  override def supportedFileAttributeViews(): util.Set[String] = HashSet("basic").asJava

  override def getSeparator: String = SEPARATOR

  override def getRootDirectories: Iterable[Path] = Seq.empty[Path].toIterable.asJavaCollection

  override def newWatchService(): WatchService = throw new UnsupportedOperationException("GCS FS does not support Watch Service")

  override def getFileStores: Iterable[FileStore] = Seq.empty[FileStore].toIterable.asJavaCollection

  override def isReadOnly: Boolean = false

  override def provider(): FileSystemProvider = gcsFileSystemProvider

  override def getPath(first: String, more: String*): Path = {
    first match {
      case gsUriRegex(chunks) => new NioGcsPath(chunks.split(SEPARATOR) ++ more.toArray[String], true)(this)
      case _ => throw new InvalidPathException(first, s"Path does not start with $PROTOCOL")
    }
  }

  /**
    * Allow instantiation of relative gcs path.
    */
  def getFlexiblePath(first: String, more: String*): Path = {
    def relativePath: Path = new NioGcsPath(first.split(SEPARATOR) ++ more.toArray[String], false)(this)
    Try(getPath(first, more: _*)) getOrElse relativePath
  }

  override def isOpen: Boolean = true

  override def close(): Unit = throw new UnsupportedOperationException("GCS FS cannot be closed")

  override def getPathMatcher(syntaxAndPattern: String): PathMatcher = FileSystems.getDefault.getPathMatcher(syntaxAndPattern)

  override def getUserPrincipalLookupService: UserPrincipalLookupService = ???
}
