package cromwell.engine.io.gcs

import java.io.File
import java.net.URI
import java.nio.file.WatchEvent.{Kind, Modifier}
import java.nio.file._
import java.util

import scala.collection.JavaConverters._
import scala.language.implicitConversions
import scala.language.postfixOps

object NioGcsPath {
  def apply(path: String)(implicit gcsFileSystem: GcsFileSystem)= {
    gcsFileSystem.getPath(path)
  }

  implicit class PathEnhanced(val path: Path) extends AnyVal {
    def asGcsPath(implicit gcsFileSystem: GcsFileSystem) = path match {
      case gcsPath: NioGcsPath => gcsPath
      case otherPath: Path => gcsFileSystem.getFlexiblePath(otherPath.toString).asInstanceOf[NioGcsPath]
      case _ => throw new IllegalArgumentException("Only GcsPaths are supported.")
    }
  }
}

class NioGcsPath(private val chunks: Array[String], absolute: Boolean)(implicit gcsFileSystem: GcsFileSystem) extends Path {
  import NioGcsPath._

  private val SEPARATOR = gcsFileSystem.getSeparator
  private val PROTOCOL = gcsFileSystem.PROTOCOL

  private val objectChunks = if(isAbsolute) chunks.tail else chunks
  private val fullPath = chunksToString(chunks)

  // Attributes
  lazy val bucket = if(isAbsolute) chunks.head else throw new UnsupportedOperationException("This Gcs path is relative, its bucket is unknown")
  val objectName = chunksToString(objectChunks)

  private def chunksToString(chunksArray: Array[String]): String = chunksArray.mkString(SEPARATOR)

  override def subpath(beginIndex: Int, endIndex: Int): Path = {
    new NioGcsPath(chunks.slice(beginIndex, endIndex), isAbsolute && beginIndex == 0)
  }

  override def toFile: File = ???

  override def resolveSibling(other: Path): Path = new NioGcsPath(chunks.init ++ other.asGcsPath.chunks, isAbsolute)

  override def resolveSibling(other: String): Path = new NioGcsPath(chunks.init ++ NioGcsPath(other).asGcsPath.chunks, isAbsolute)

  override def getFileSystem: FileSystem = gcsFileSystem

  override def getName(index: Int): Path = new NioGcsPath(Array(chunks(index)), isAbsolute && index == 0)

  override def getParent: Path = new NioGcsPath(chunks.init, isAbsolute)

  override def toAbsolutePath: Path = if (isAbsolute) this else throw new UnsupportedOperationException("Cannot make a relative GCS path absolute.")

  override def relativize(other: Path): Path = ???

  override def getNameCount: Int = chunks.length

  override def toUri: URI = new URI(toString)

  override def compareTo(other: Path): Int = ???

  override def register(watcher: WatchService, events: Array[Kind[_]], modifiers: Modifier*): WatchKey = ???

  override def register(watcher: WatchService, events: Kind[_]*): WatchKey = ???

  override def getFileName: Path = new NioGcsPath(Array(chunks.last), false)

  override def getRoot: Path = new NioGcsPath(Array(bucket), true)

  override def iterator(): util.Iterator[Path] = (chunks map { gcsFileSystem.getFlexiblePath(_) } iterator).asJava

  override def normalize(): Path = if (isAbsolute) this else throw new UnsupportedOperationException("Cannot normalize a relative GCS path.")

  override def endsWith(other: Path): Boolean = chunks.endsWith(other.asGcsPath.chunks)

  override def endsWith(other: String): Boolean = chunks.endsWith(gcsFileSystem.getFlexiblePath(other).asGcsPath.chunks)

  override def resolve(other: Path): Path = {
    if (other.isAbsolute) other
    else new NioGcsPath(chunks ++ other.asGcsPath.chunks, isAbsolute)
  }

  override def resolve(other: String): Path = {
    val otherPath = gcsFileSystem.getFlexiblePath(other)
    if (otherPath.isAbsolute) otherPath
    else new NioGcsPath(chunks ++ otherPath.asGcsPath.chunks, isAbsolute)
  }

  override def toRealPath(options: LinkOption*): Path = this

  override def startsWith(other: Path): Boolean = chunks.startsWith(other.asGcsPath.chunks)

  override def startsWith(other: String): Boolean = chunks.startsWith(gcsFileSystem.getFlexiblePath(other).asGcsPath.chunks)

  override def toString: String = s"$PROTOCOL$fullPath"

  override def isAbsolute: Boolean = absolute
}
