package cromwell.core.path.proxy

import java.io.File
import java.net.URI
import java.nio.file.WatchEvent.{Kind, Modifier}
import java.nio.file._
import java.util

import scala.util.Try

class PathProxy(delegate: Path, injectedFileSystem: FileSystem) extends Path {
  def unbox[T](clazz: Class[T]): Try[T] = Try {
    clazz.cast(delegate)
  }

  override def getFileSystem: FileSystem = injectedFileSystem

  /* delegated */
  override def subpath(beginIndex: Int, endIndex: Int): Path = delegate.subpath(beginIndex, endIndex)
  override def toFile: File = delegate.toFile
  override def resolveSibling(other: Path): Path = delegate.resolveSibling(other)
  override def resolveSibling(other: String): Path = delegate.resolveSibling(other)
  override def isAbsolute: Boolean = delegate.isAbsolute
  override def getName(index: Int): Path = delegate.getName(index)
  override def getParent: Path = delegate.getParent
  override def toAbsolutePath: Path = delegate.toAbsolutePath
  override def relativize(other: Path): Path = delegate.relativize(other)
  override def getNameCount: Int = delegate.getNameCount
  override def toUri: URI = delegate.toUri
  override def compareTo(other: Path): Int = delegate.compareTo(other)
  override def register(watcher: WatchService, events: Array[Kind[_]], modifiers: Modifier*): WatchKey = delegate.register(watcher, events, modifiers: _*)
  override def register(watcher: WatchService, events: Kind[_]*): WatchKey = delegate.register(watcher, events: _*)
  override def getFileName: Path = delegate.getFileName
  override def getRoot: Path = delegate.getRoot
  override def iterator(): util.Iterator[Path] = delegate.iterator()
  override def normalize(): Path = delegate.normalize()
  override def endsWith(other: Path): Boolean = delegate.endsWith(other)
  override def endsWith(other: String): Boolean = delegate.endsWith(other)
  override def resolve(other: Path): Path = delegate.resolve(other)
  override def resolve(other: String): Path = delegate.resolve(other)
  override def startsWith(other: Path): Boolean = delegate.startsWith(other)
  override def startsWith(other: String): Boolean = delegate.startsWith(other)
  override def toRealPath(options: LinkOption*): Path = delegate.toRealPath(options: _*)
}
