package cromwell.engine.io.gcs

import java.lang.Iterable
import java.nio.file._
import java.nio.file.attribute.UserPrincipalLookupService
import java.nio.file.spi.FileSystemProvider
import java.util.{Collections, Set => JSet}

import cromwell.engine.PathString

import scala.concurrent.ExecutionContext

object GcsFileSystem {
  import PathString._
  def instance(interface: GoogleCloudStorage, root: String)(implicit executionContext: ExecutionContext) = {
    if (root.isGcsUrl) new GcsFileSystem(GcsFileSystemProvider.instance(interface), root)
    else throw new IllegalArgumentException(s"$root is not am absolute GCS path")
  }
  val Separator = "/"
  private [io] val Protocol = "gs://"
  private val GsUriRegex = s"""$Protocol(.*)""".r
  private val AttributeViews = Collections.singleton("basic")
}

/**
  * Implements the java.nio.FileSystem interface for GoogleCloudStorage.
  *
  */
class GcsFileSystem private (gcsFileSystemProvider: GcsFileSystemProvider, gcsRoot: String) extends FileSystem {

  import GcsFileSystem._

  val root = gcsRoot match {
    case GsUriRegex(chunks) => new NioGcsPath(chunks.split(Separator), true)(this)
    case _ => throw new InvalidPathException(gcsRoot, s"Root of GCS file system must be an absolute GCS path: $gcsRoot")
  }

  override def supportedFileAttributeViews(): JSet[String] = AttributeViews
  override def getSeparator: String = Separator
  override def getRootDirectories: Iterable[Path] = Collections.singleton(root)
  override def newWatchService(): WatchService = throw new NotImplementedError("GCS FS does not support Watch Service at this time")
  override def getFileStores: Iterable[FileStore] = Collections.emptyList()
  override def isReadOnly: Boolean = false
  override def provider(): FileSystemProvider = gcsFileSystemProvider
  override def isOpen: Boolean = true
  override def close(): Unit = throw new UnsupportedOperationException("GCS FS cannot be closed")
  override def getPathMatcher(syntaxAndPattern: String): PathMatcher = FileSystems.getDefault.getPathMatcher(syntaxAndPattern)
  override def getUserPrincipalLookupService: UserPrincipalLookupService = throw new UnsupportedOperationException()
  override def getPath(first: String, more: String*): Path = {
    first match {
      case GsUriRegex(chunks) =>
        val path = new NioGcsPath(chunks.split(Separator) ++ more.toArray[String], true)(this)
        if (path.startsWith(root)) path
        else throw new InvalidPathException(first, s"Path $path has a different root than this filesystem: $root")
      case empty if empty.isEmpty => new NioGcsPath(Array.empty[String] ++ more.toArray[String], false)(this)
      case _ => new NioGcsPath(first.split(Separator) ++ more.toArray[String], false)(this)
    }
  }
}
