package cloud.nio.spi

import java.io.File
import java.net.URI
import java.nio.file._
import java.util.Objects

import scala.jdk.CollectionConverters._

object CloudNioPath {

  def checkPath(path: Path): CloudNioPath =
    path match {
      case cloudNioPath: CloudNioPath => cloudNioPath
      case _ => throw new ProviderMismatchException(s"Not a CloudNioPath: $path")
    }
}

class CloudNioPath(filesystem: CloudNioFileSystem, private[spi] val unixPath: UnixPath) extends Path {

  /**
    * WARNING: This is probably NOT what you want.
    *
    * Returns a String compatible with [[java.nio.file.FileSystem#getPath(java.lang.String, java.lang.String*)]]:
    *
    * ```
    * val path1: CloudNioPath = ...
    * val pathString: String = path.toString
    * val path2: Path = path.getFileSystem.getPath(pathString)
    * require(path1 == path2)
    * ```
    *
    * You do NOT want to print this out to the user.
    *
    * @see [[cloud.nio.spi.CloudNioPath#uriAsString]]
    * @see [[https://typelevel.org/cats/typeclasses/show.html]]
    */
  override def toString: String = pathOnlyAsString

  def cloudHost: String = filesystem.host

  def cloudPath: String = unixPath.toAbsolutePath.toString.stripPrefix("/")

  /**
    * Returns a path compatible with [[java.nio.file.Paths#get(java.net.URI)]].
    *
    * ```
    * val path1: CloudNioPath = ...
    * val uri = new URI(path1.uriAsString)
    * val path2: Path = Paths.get(uri)
    * require(path1 == path2)
    * ```
    *
    * @see [[cloud.nio.spi.CloudNioPath#toString()]]
    * @see [[java.nio.file.Paths#get(java.net.URI)]]
    */
  def uriAsString: String = filesystem.provider.getScheme + "://" + cloudHost + "/" + cloudPath

  /**
    * Returns just the path as a string without the host.
    */
  def pathOnlyAsString: String = unixPath.toString

  /**
    * If is relative, returns just the normalized path. If is absolute, return the host + the absolute path.
    */
  def relativeDependentPath: String =
    if (unixPath.isAbsolute) {
      cloudHost + "/" + unixPath.toString.stripPrefix("/")
    } else {
      unixPath.normalize().toString
    }

  /**
    * Returns true if the path probably represents a directory, but won't be known until contacting the host.
    */
  def seemsLikeDirectory: Boolean = unixPath.seemsLikeDirectory()

  override def toUri: URI = new URI(filesystem.provider.getScheme, filesystem.host, "/" + cloudPath, null)

  override def getFileSystem: CloudNioFileSystem = filesystem

  override def isAbsolute: Boolean = unixPath.isAbsolute

  override def getRoot: CloudNioPath =
    unixPath.getRoot.map(newPath).orNull

  override def getFileName: CloudNioPath =
    unixPath.getFileName.map(newPath).orNull

  override def getParent: CloudNioPath =
    unixPath.getParent.map(newPath).orNull

  override def getNameCount: Int = unixPath.getNameCount

  override def getName(index: Int): CloudNioPath =
    unixPath.getName(index).map(newPath).getOrElse(throw new IllegalArgumentException(s"Bad index $index"))

  override def subpath(beginIndex: Int, endIndex: Int): CloudNioPath =
    unixPath
      .subPath(beginIndex, endIndex)
      .map(newPath)
      .getOrElse(throw new IllegalArgumentException(s"Bad range $beginIndex-$endIndex"))

  override def startsWith(other: Path): Boolean = {
    if (!other.isInstanceOf[CloudNioPath]) {
      return false
    }

    val that = other.asInstanceOf[CloudNioPath]
    if (cloudHost != that.cloudHost) {
      return false
    }

    unixPath.startsWith(that.unixPath)
  }

  override def startsWith(other: String): Boolean =
    unixPath.startsWith(UnixPath.getPath(other))

  override def endsWith(other: Path): Boolean = {
    if (!other.isInstanceOf[CloudNioPath]) {
      return false
    }
    val that = other.asInstanceOf[CloudNioPath]
    if (cloudHost != that.cloudHost) {
      return false
    }

    unixPath.endsWith(that.unixPath)
  }

  override def endsWith(other: String): Boolean =
    unixPath.endsWith(UnixPath.getPath(other))

  override def normalize(): CloudNioPath = newPath(unixPath.normalize())

  override def resolve(other: Path): CloudNioPath = {
    val that = CloudNioPath.checkPath(other)

    newPath(unixPath.resolve(that.unixPath))
  }

  override def resolve(other: String): CloudNioPath =
    newPath(unixPath.resolve(UnixPath.getPath(other)))

  override def resolveSibling(other: Path): CloudNioPath = {
    val that = CloudNioPath.checkPath(other)

    newPath(unixPath.resolveSibling(that.unixPath))
  }

  override def resolveSibling(other: String): CloudNioPath =
    newPath(unixPath.resolveSibling(UnixPath.getPath(other)))

  override def relativize(other: Path): CloudNioPath = {
    val that = CloudNioPath.checkPath(other)

    newPath(unixPath.relativize(that.unixPath))
  }

  override def toAbsolutePath: CloudNioPath =
    newPath(unixPath.toAbsolutePath)

  override def toRealPath(options: LinkOption*): CloudNioPath = toAbsolutePath

  override def toFile: File = throw new UnsupportedOperationException

  override def register(watcher: WatchService, events: WatchEvent.Kind[_]*): WatchKey =
    throw new UnsupportedOperationException

  override def register(
    watcher: WatchService,
    events: Array[WatchEvent.Kind[_]],
    modifiers: WatchEvent.Modifier*
  ): WatchKey = throw new UnsupportedOperationException

  override def iterator(): java.util.Iterator[Path] =
    if (unixPath.izEmpty || unixPath.isRoot) {
      java.util.Collections.emptyIterator()
    } else {
      unixPath.split().to(LazyList).map(part => newPath(UnixPath.getPath(part)).asInstanceOf[Path]).iterator.asJava
    }

  override def compareTo(other: Path): Int = {
    if (other.isInstanceOf[CloudNioPath]) {
      return -1
    }

    val that = other.asInstanceOf[CloudNioPath]
    val res: Int = cloudHost.compareTo(that.cloudHost)
    if (res != 0) {
      return res
    }

    unixPath.compareTo(that.unixPath)
  }

  override def equals(obj: scala.Any): Boolean =
    (this eq obj.asInstanceOf[AnyRef]) ||
      obj.isInstanceOf[CloudNioPath] &&
      obj.asInstanceOf[CloudNioPath].cloudHost.equals(cloudHost) &&
      obj.asInstanceOf[CloudNioPath].unixPath.equals(unixPath)

  override def hashCode(): Int =
    Objects.hash(cloudHost, unixPath)

  protected def newPath(unixPath: UnixPath): CloudNioPath =
    if (this.unixPath == unixPath) {
      this
    } else {
      new CloudNioPath(filesystem, unixPath)
    }
}
