package cromwell.filesystems.gcs

import java.lang.Iterable
import java.nio.file._
import java.nio.file.attribute.UserPrincipalLookupService
import java.nio.file.spi.FileSystemProvider
import java.util.{Collections, Set => JSet}

import scala.language.postfixOps

case class NotAGcsPathException(path: String) extends IllegalArgumentException(s"$path is not a valid GCS path.")

object GcsFileSystem {
  val Separator = "/"
  private[gcs] val Protocol = "gs://"
  private val GsUriRegex = s"""$Protocol(.*)""".r
  private val AttributeViews = Collections.singleton("basic")

  val defaultGcsFileSystem = GcsFileSystemProvider.defaultProvider.getFileSystem

  def isAbsoluteGcsPath(str: String) = str match {
    case GsUriRegex(chunks) => true
    case _ => false
  }

  def apply(provider: GcsFileSystemProvider) = new GcsFileSystem(provider)
}

/**
  * Implements the java.nio.FileSystem interface for GoogleCloudStorage.
  */
class GcsFileSystem private(val gcsFileSystemProvider: GcsFileSystemProvider) extends FileSystem {

  import GcsFileSystem._

  override def supportedFileAttributeViews(): JSet[String] = AttributeViews

  override def getSeparator: String = Separator

  override def getRootDirectories: Iterable[Path] = Collections.emptyList[Path]

  override def newWatchService(): WatchService = throw new NotImplementedError("GCS FS does not support Watch Service at this time")

  override def getFileStores: Iterable[FileStore] = Collections.emptyList()

  override def isReadOnly: Boolean = false

  override def provider(): FileSystemProvider = gcsFileSystemProvider

  override def isOpen: Boolean = true

  override def close(): Unit = throw new UnsupportedOperationException("GCS FS cannot be closed")

  override def getPathMatcher(syntaxAndPattern: String): PathMatcher = FileSystems.getDefault.getPathMatcher(syntaxAndPattern)

  override def getUserPrincipalLookupService: UserPrincipalLookupService = throw new UnsupportedOperationException()

  private def buildPath(first: String, more: Seq[String], forceDirectory: Boolean) = {
    val directory = forceDirectory || (more.isEmpty && first.endsWith(Separator)) || more.lastOption.exists(_.endsWith(Separator))
    first match {
      case GsUriRegex(chunks) => new NioGcsPath(chunks.split(Separator) ++ more.toArray[String], true, directory)(this)
      case empty if empty.isEmpty => new NioGcsPath(Array.empty[String] ++ more.toArray[String], false, false)(this)
      case _ => throw new NotAGcsPathException(s"$first is not a gcs path")
    }
  }

  override def getPath(first: String, more: String*): Path = buildPath(first, more, forceDirectory = false)

  def getPathAsDirectory(first: String, more: String*): Path = buildPath(first, more, forceDirectory = true)
}
