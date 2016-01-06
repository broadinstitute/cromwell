package cromwell.engine.io.gcs

import java.lang.Iterable
import java.nio.file._
import java.nio.file.attribute.UserPrincipalLookupService
import java.nio.file.spi.FileSystemProvider
import java.util.{Collections, Set => JSet}

import cromwell.engine.PathString

import scala.concurrent.ExecutionContext
import scala.util.Try

object GcsFileSystem {
  import PathString._
  def getInstance(interface: GoogleCloudStorage, root: String)(implicit executionContext: ExecutionContext) = {
    if(root.isGcsUrl) new GcsFileSystem(GcsFileSystemProvider.getInstance(interface), root)
    else throw new IllegalArgumentException(s"$root is not am absolute GCS path")
  }
  private val SEPARATOR = "/"
  private [io] val PROTOCOL = "gs://"
  private val GS_URI_REGEX = s"""$PROTOCOL(.*)""".r
  private val attributeViews = Collections.singleton("basic")
}

/**
  * Implements the java.nio.FileSystem interface for GoogleCloudStorage.
  */
class GcsFileSystem private (gcsFileSystemProvider: GcsFileSystemProvider, gcsRoot: String) extends FileSystem {

  import GcsFileSystem._

  val root = getPath(gcsRoot)

  override def supportedFileAttributeViews(): JSet[String] = attributeViews

  override def getSeparator: String = SEPARATOR

  override def getRootDirectories: Iterable[Path] = Collections.singleton(root)

  override def newWatchService(): WatchService = throw new UnsupportedOperationException("GCS FS does not support Watch Service")

  override def getFileStores: Iterable[FileStore] = Collections.emptyList()

  override def isReadOnly: Boolean = false

  override def provider(): FileSystemProvider = gcsFileSystemProvider

  override def getPath(first: String, more: String*): Path = {
    first match {
      case GS_URI_REGEX(chunks) => new NioGcsPath(chunks.split(SEPARATOR) ++ more.toArray[String], true)(this)
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

  override def getUserPrincipalLookupService: UserPrincipalLookupService = throw new NotImplementedError()
}
