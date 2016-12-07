package cromwell.core.path.proxy

import java.net.URI
import java.nio.channels.SeekableByteChannel
import java.nio.file.DirectoryStream.Filter
import java.nio.file._
import java.nio.file.attribute.{BasicFileAttributes, FileAttribute, FileAttributeView}
import java.nio.file.spi.FileSystemProvider
import java.util

import akka.actor.ActorSystem
import cromwell.core.path.CustomRetryParams
import cromwell.core.retry.Retry

import scala.concurrent.{Await, Future}

class RetryableFileSystemProviderProxy[T <: FileSystemProvider](delegate: T, retryParams: CustomRetryParams = CustomRetryParams.Default)(implicit actorSystem: ActorSystem) extends FileSystemProvider {
  private val iOExecutionContext = actorSystem.dispatchers.lookup("akka.dispatchers.io-dispatcher")

  // the nio interface is synchronous so we need to wait for the result
  def withRetry[U](f: () => U): U = Await.result(
    Retry.withRetry(
      () => Future(f())(iOExecutionContext),
      retryParams.maxRetries,
      retryParams.backoff,
      retryParams.isTransient,
      retryParams.isFatal
    ),
    retryParams.timeout
  )

  override def getPath(uri: URI): Path = {
    val path = delegate.getPath(uri)
    new PathProxy(path, new FileSystemProxy(path.getFileSystem, this))
  }
  override def newFileSystem(uri: URI, env: util.Map[String, _]): FileSystem = {
    new FileSystemProxy(delegate.newFileSystem(uri, env), this)
  }
  override def getScheme: String = delegate.getScheme
  override def getFileSystem(uri: URI): FileSystem = {
    new FileSystemProxy(delegate.getFileSystem(uri), this)
  }
  override def getFileStore(path: Path): FileStore = delegate.getFileStore(path)

  /* retried operations */
  override def move(source: Path, target: Path, options: CopyOption*): Unit = withRetry { () => delegate.move(source, target, options: _*) }
  override def checkAccess(path: Path, modes: AccessMode*): Unit = withRetry { () => delegate.checkAccess(path, modes: _*) }
  override def createDirectory(dir: Path, attrs: FileAttribute[_]*): Unit = withRetry { () => delegate.createDirectory(dir, attrs: _*) }
  override def newByteChannel(path: Path, options: util.Set[_ <: OpenOption], attrs: FileAttribute[_]*): SeekableByteChannel = withRetry { () => delegate.newByteChannel(path, options, attrs: _*) }
  override def isHidden(path: Path): Boolean = withRetry { () => delegate.isHidden(path) }
  override def copy(source: Path, target: Path, options: CopyOption*): Unit = withRetry { () => delegate.copy(source, target, options: _*) }
  override def delete(path: Path): Unit = withRetry { () => delegate.delete(path) }
  override def newDirectoryStream(dir: Path, filter: Filter[_ >: Path]): DirectoryStream[Path] = withRetry { () => delegate.newDirectoryStream(dir, filter) }
  override def setAttribute(path: Path, attribute: String, value: scala.Any, options: LinkOption*): Unit = withRetry { () => delegate.setAttribute(path, attribute, value, options: _*) }
  override def readAttributes[A <: BasicFileAttributes](path: Path, `type`: Class[A], options: LinkOption*): A = withRetry { () => delegate.readAttributes(path, `type`, options: _*) }
  override def readAttributes(path: Path, attributes: String, options: LinkOption*): util.Map[String, AnyRef] = withRetry { () => delegate.readAttributes(path, attributes, options: _*) }
  override def isSameFile(path: Path, path2: Path): Boolean = withRetry { () => delegate.isSameFile(path, path2) }
  override def getFileAttributeView[V <: FileAttributeView](path: Path, `type`: Class[V], options: LinkOption*): V = withRetry { () => delegate.getFileAttributeView(path, `type`, options: _*) }
}
