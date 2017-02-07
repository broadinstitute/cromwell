package cromwell.core.path

import java.net.URI
import java.nio.file.WatchEvent.{Kind, Modifier}
import java.nio.file.{LinkOption, WatchKey, WatchService}

import scala.collection.JavaConverters._

/**
  * Implements methods with the same names and signatures as java.nio.Path
  *
  * @see [[cromwell.core.path.Path]]
  * @see [[cromwell.core.path.BetterFileMethods]]
  * @see [[cromwell.core.path.EvenBetterPathMethods]]
  */
trait NioPathMethods {
  self: Path =>

  final def subpath(beginIndex: Int, endIndex: Int): Path = newPath(nioPathPrivate.subpath(beginIndex, endIndex))

  final def toFile: java.io.File = nioPathPrivate.toFile

  final def resolveSibling(other: Path): Path = newPath(nioPathPrivate.resolveSibling(other.nioPathPrivate))

  final def resolveSibling(other: String): Path = newPath(nioPathPrivate.resolveSibling(other))

  final def isAbsolute: Boolean = nioPathPrivate.isAbsolute

  final def getName(index: Int): Path = newPath(nioPathPrivate.getName(index))

  final def getParent: Path = newPathOrNull(nioPathPrivate.getParent)

  final def toAbsolutePath: Path = newPath(nioPathPrivate.toAbsolutePath)

  final def relativize(other: Path): Path = newPath(nioPathPrivate.relativize(other.nioPathPrivate))

  final def getNameCount: Int = nioPathPrivate.getNameCount

  final def toUri: URI = nioPathPrivate.toUri

  final def compareTo(other: Path): Int = nioPathPrivate.compareTo(other.nioPathPrivate)

  final def register(watcher: WatchService, events: Array[Kind[_]], modifiers: Modifier*): WatchKey =
    nioPathPrivate.register(watcher, events, modifiers: _*)

  final def register(watcher: WatchService, events: Kind[_]*): WatchKey = nioPathPrivate.register(watcher, events: _*)

  final def getFileName: Path = newPathOrNull(nioPathPrivate.getFileName)

  final def getRoot: Path = newPathOrNull(nioPathPrivate.getRoot)

  final def iterator(): java.util.Iterator[Path] = nioPathPrivate.iterator().asScala.map(newPath).asJava

  final def normalize(): Path = newPath(nioPathPrivate.normalize())

  final def endsWith(other: Path): Boolean = nioPathPrivate.endsWith(other.nioPathPrivate)

  final def endsWith(other: String): Boolean = nioPathPrivate.endsWith(other)

  final def resolve(other: Path): Path = newPath(nioPathPrivate.resolve(other.nioPathPrivate))

  final def resolve(other: String): Path = newPath(nioPathPrivate.resolve(other))

  final def startsWith(other: Path): Boolean = nioPathPrivate.startsWith(other.nioPathPrivate)

  final def startsWith(other: String): Boolean = nioPathPrivate.startsWith(other)

  final def toRealPath(options: LinkOption*): Path = newPath(nioPathPrivate.toRealPath(options: _*))
}
