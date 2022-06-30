package cromwell.filesystems.oss.nio

import java.io.File
import java.net.URI
import java.nio.file._
import java.util
import java.util.Objects

import com.google.common.collect.UnmodifiableIterator

object OssStoragePath {
  def checkOssStoragePath(other: Path): OssStoragePath = {
    if (!other.isInstanceOf[OssStoragePath]) {
      throw new ProviderMismatchException(s"Not a oss storage path $other")
    }

    other.asInstanceOf[OssStoragePath]
  }

  def getPath(filesystem: OssStorageFileSystem, path: UnixPath) = new OssStoragePathImpl(filesystem, path)

  def getPath(filesystem: OssStorageFileSystem, first: String, more: String*) = new OssStoragePathImpl(filesystem, UnixPath.getPath(first, more: _*))

  def checkPath(path: Path): OssStoragePath = {
    if (!path.isInstanceOf[OssStoragePath]) {
      throw new ProviderMismatchException(s"Not an oss storage path $path")
    }

    path.asInstanceOf[OssStoragePath]
  }

}

trait OssStoragePath extends Path {
  def bucket = ""

  def key = ""

  def path: UnixPath = UnixPath.EMPTY_PATH

  def seemsLikeDirectory = false

  def pathAsString: String = ""

  override def getFileSystem: OssStorageFileSystem = throw new UnsupportedOperationException

  override def isAbsolute: Boolean = throw new UnsupportedOperationException

  override def getRoot: OssStoragePath = throw new UnsupportedOperationException

  override def getFileName: OssStoragePath = throw new UnsupportedOperationException

  override def getParent: OssStoragePath = throw new UnsupportedOperationException

  override def getNameCount: Int = throw new UnsupportedOperationException

  override def getName(index: Int): OssStoragePath = throw new UnsupportedOperationException

  override def subpath(beginIndex: Int, endIndex: Int): OssStoragePath = throw new UnsupportedOperationException

  override def startsWith(other: Path): Boolean = throw new UnsupportedOperationException

  override def startsWith(other: String): Boolean = throw new UnsupportedOperationException

  override def endsWith(other: Path): Boolean = throw new UnsupportedOperationException

  override def endsWith(other: String): Boolean = throw new UnsupportedOperationException

  override def normalize(): OssStoragePath = throw new UnsupportedOperationException

  override def resolve(other: Path): OssStoragePath = throw new UnsupportedOperationException

  override def resolve(other: String): OssStoragePath = throw new UnsupportedOperationException

  override def resolveSibling(other: Path): OssStoragePath = throw new UnsupportedOperationException

  override def resolveSibling(other: String): OssStoragePath = throw new UnsupportedOperationException

  override def relativize(other: Path): OssStoragePath = throw new UnsupportedOperationException

  override def toAbsolutePath: OssStoragePath = throw new UnsupportedOperationException

  override def toRealPath(options: LinkOption*): OssStoragePath = throw new UnsupportedOperationException

  override def toFile: File = throw new UnsupportedOperationException

  override def register(watcher: WatchService, events: WatchEvent.Kind[_]*): WatchKey = throw new UnsupportedOperationException

  override def register(watcher: WatchService, events: Array[WatchEvent.Kind[_]], modifiers: WatchEvent.Modifier*): WatchKey = throw new UnsupportedOperationException

  override def iterator(): util.Iterator[Path] = throw new UnsupportedOperationException

  override def compareTo(other: Path): Int = throw new UnsupportedOperationException

  override def toUri: URI = throw new UnsupportedOperationException
}

final case class OssStoragePathImpl(filesystem: OssStorageFileSystem, override val path: UnixPath = UnixPath.EMPTY_PATH) extends  OssStoragePath {

  override def pathAsString: String = toUri.toString

  override def bucket: String = filesystem.bucket

  override def key: String = toAbsolutePath.toString.stripPrefix("/")

  override def getFileSystem: OssStorageFileSystem = filesystem

  override def isAbsolute: Boolean = path.isAbsolute

  override def getRoot: OssStoragePath = path.getRoot map {path => newPath(path)} getOrElse NullOssStoragePath(filesystem)

  override def getFileName: OssStoragePath = path.getFileName map {path => newPath(path)} getOrElse NullOssStoragePath(filesystem)

  override def getParent: OssStoragePath = path.getParent map {path => newPath(path)} getOrElse NullOssStoragePath(filesystem)

  override def getNameCount: Int = path.getNameCount

  override def getName(index: Int): OssStoragePath = path.getName(index) map {path => newPath(path)} getOrElse NullOssStoragePath(filesystem)

  override def subpath(beginIndex: Int, endIndex: Int): OssStoragePath = path.subPath(beginIndex, endIndex) map {path => newPath(path)} getOrElse NullOssStoragePath(filesystem)

  override def startsWith(other: Path): Boolean = {
    if (!other.isInstanceOf[OssStoragePath]) {
      return false
    }

    val that = other.asInstanceOf[OssStoragePath]
    if (bucket != that.bucket) {
      return false
    }

    path.startsWith(that.path)
  }

  override def startsWith(other: String): Boolean = {
    path.startsWith(UnixPath.getPath(other))
  }

  override def endsWith(other: Path): Boolean = {
    if (!other.isInstanceOf[OssStoragePath]) {
      return false
    }
    val that = other.asInstanceOf[OssStoragePath]
    if (bucket != that.bucket) {
      return false
    }

    path.endsWith(that.path)
  }

  override def endsWith(other: String): Boolean = {
    path.endsWith(UnixPath.getPath(other))
  }

  override def normalize(): OssStoragePath = newPath(path.normalize())

  override def resolve(other: Path): OssStoragePath = {
    val that = OssStoragePath.checkOssStoragePath(other)

    newPath(path.resolve(that.path))
  }

  override def resolve(other: String): OssStoragePath = {
    newPath(path.resolve(UnixPath.getPath(other)))
  }

  override def resolveSibling(other: Path): OssStoragePath = {
    val that = OssStoragePath.checkOssStoragePath(other)

    newPath(path.resolveSibling(that.path))
  }

  override def resolveSibling(other: String): OssStoragePath = {
    newPath(path.resolveSibling(UnixPath.getPath(other)))
  }

  override def relativize(other: Path): OssStoragePath = {
    val that = OssStoragePath.checkOssStoragePath(other)

    newPath(path.relativize(that.path))
  }

  /**
   * currently a mocked one
   */
  override def toAbsolutePath: OssStoragePath = {
    newPath(path.toAbsolutePath())
  }

  override def toRealPath(options: LinkOption*): OssStoragePath = toAbsolutePath

  override def toFile: File = throw new UnsupportedOperationException

  override def register(watcher: WatchService, events: WatchEvent.Kind[_]*): WatchKey = throw new UnsupportedOperationException

  override def register(watcher: WatchService, events: Array[WatchEvent.Kind[_]], modifiers: WatchEvent.Modifier*): WatchKey = throw new UnsupportedOperationException

  override def iterator(): util.Iterator[Path] = {
    if (path.isEmpty() || path.isRoot) {
      return util.Collections.emptyIterator()
    }

    new PathIterator()
  }

  override def compareTo(other: Path): Int = {
    if (other.isInstanceOf[OssStoragePath]) {
      return -1
    }

    val that = other.asInstanceOf[OssStoragePath]
    val res: Int = bucket.compareTo(that.bucket)
    if (res != 0) {
      return res
    }

    path.compareTo(that.path)
  }

  override def seemsLikeDirectory = path.seemsLikeDirectory()

  override def equals(obj: scala.Any): Boolean = {
    (this eq obj.asInstanceOf[AnyRef]) || obj.isInstanceOf[OssStoragePath] && obj.asInstanceOf[OssStoragePath].bucket.equals(bucket) && obj.asInstanceOf[OssStoragePath].path.equals(path)
  }

  override def hashCode(): Int = {
    Objects.hash(bucket, toAbsolutePath.toString)
  }

  override def toString: String = path.toString

  override def toUri: URI = new URI("oss", bucket, toAbsolutePath.toString, None.orNull)

  private[this] def newPath(unixPath: UnixPath): OssStoragePath = {
    if (unixPath == path) {
      this
    } else {
      OssStoragePathImpl(filesystem, unixPath)
    }
  }


  class PathIterator extends UnmodifiableIterator[Path] {
    val delegate = path.split()

    override def next(): OssStoragePath = newPath(UnixPath.getPath(delegate.next()))

    override def hasNext: Boolean = delegate.hasNext
  }
}

final case class NullOssStoragePath(filesystem: OssStorageFileSystem) extends OssStoragePath
