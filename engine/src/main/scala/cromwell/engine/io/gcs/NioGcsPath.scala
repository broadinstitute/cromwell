package cromwell.engine.io.gcs

import java.io.File
import java.net.URI
import java.nio.file.WatchEvent.{Kind, Modifier}
import java.nio.file._
import java.util

import scala.collection.JavaConverters._
import scala.language.implicitConversions
import scala.language.postfixOps
import scala.util.Try

object NioGcsPath {
  def apply(path: String)(implicit gcsFileSystem: GcsFileSystem)= gcsFileSystem.getPath(path)

  implicit class PathEnhanced(val path: Path) extends AnyVal {
    def asGcsPath(implicit gcsFileSystem: GcsFileSystem) = path match {
      case gcsPath: NioGcsPath => gcsPath
      case otherPath: Path => gcsFileSystem.getPath(otherPath.toString).asInstanceOf[NioGcsPath]
      case _ => throw new IllegalArgumentException("Only GcsPaths are supported.")
    }
  }

  val Protocol = GcsFileSystem.Protocol
}

/**
  * NOTE: Currently called NioGcsPath so it can exist alongside the current GcsPath class.
  *       If this approach was to be validated the current GcsPath class would be replaced by this one.
  * This class proposes an implementation of the java.nio.Path interface for GoogleCloudStorage.
  * The following methods are yet to be implemented:
  *   relativize
  *   compareTo
  * @param chunks array containing all parts of the path in between separators - except the protocol (gs://)
  *               eg: gs://path/to/resource.txt -> chunks = [path, to, resource.txt]
  * @param absolute true if this path is to be considered absolute.
  *                 Only absolute GCS paths can be used to actually locate resources.
  *                 Calling methods on an absolute path can return a relative paths (eg subpath).
  * @param gcsFileSystem the gcsFileSystem to be used when performing operations on this path
  */
class NioGcsPath(private val chunks: Array[String], absolute: Boolean)(implicit gcsFileSystem: GcsFileSystem) extends Path {
  import NioGcsPath._

  private val separator = GcsFileSystem.Separator

  private val objectChunks = chunks match {
    case values if isAbsolute && values.nonEmpty => values.tail
    case _ => chunks
  }

  private val fullPath = chunksToString(chunks)

  lazy val bucket: String = chunks match {
    case values if values.isEmpty && isAbsolute => throw new IllegalStateException("An absolute gcs path cannot be empty")
    case _ => if(isAbsolute) chunks.head else gcsFileSystem.root.asGcsPath.bucket
  }

  val objectName = chunksToString(objectChunks)

  private def chunksToString(chunksArray: Array[String]): String = chunksArray.mkString(separator)

  override def subpath(beginIndex: Int, endIndex: Int): Path = {
    new NioGcsPath(chunks.slice(beginIndex, endIndex), isAbsolute && beginIndex == 0)
  }

  override def toFile: File = throw new UnsupportedOperationException("A GCS path cannot be converted to a File.")

  override def resolveSibling(other: Path): Path = new NioGcsPath(getParent.asGcsPath.chunks ++ other.asGcsPath.chunks, isAbsolute)

  override def resolveSibling(other: String): Path = new NioGcsPath(getParent.asGcsPath.chunks ++ gcsFileSystem.getPath(other).asGcsPath.chunks, isAbsolute)

  override def getFileSystem: FileSystem = gcsFileSystem

  override def getName(index: Int): Path = new NioGcsPath(Array(chunks(index)), isAbsolute && index == 0)

  override def getParent: Path = chunks match {
    case values if values.isEmpty || values.length == 1 => null
    case values => new NioGcsPath(values.init, isAbsolute)
  }

  override def toAbsolutePath: Path = if (isAbsolute) this else gcsFileSystem.root.resolve(this)

  override def relativize(other: Path): Path = throw new NotImplementedError()

  override def getNameCount: Int = chunks.length

  override def toUri: URI = throw new UnsupportedOperationException()

  override def compareTo(other: Path): Int = throw new NotImplementedError()

  override def register(watcher: WatchService, events: Array[Kind[_]], modifiers: Modifier*): WatchKey = throw new UnsupportedOperationException()

  override def register(watcher: WatchService, events: Kind[_]*): WatchKey = throw new UnsupportedOperationException()

  override def getFileName: Path = chunks match {
    case values if values.isEmpty => null
    case _ => new NioGcsPath(Array(chunks.last), isAbsolute && chunks.length == 1)
  }

  override def getRoot: Path = new NioGcsPath(Array(bucket), true)

  override def iterator(): util.Iterator[Path] = (chunks map { elt => new NioGcsPath(Array(elt), false).asInstanceOf[Path] } iterator).asJava

  override def normalize(): Path = if (isAbsolute) this else throw new UnsupportedOperationException("Cannot normalize a relative GCS path.")

  override def endsWith(other: Path): Boolean = {
    other match {
      case rel: NioGcsPath if !isAbsolute && rel.isAbsolute => false
      case _: NioGcsPath => chunks.endsWith(other.asGcsPath.chunks)
      case _ => false
    }
  }

  override def endsWith(other: String): Boolean = {
    Try(gcsFileSystem.getPath(other)) map {
      case rel: NioGcsPath if !isAbsolute && rel.isAbsolute => false
      case path@(_: NioGcsPath) => chunks.endsWith(path.asGcsPath.chunks)
      case _ => false
    } getOrElse false
  }

  override def resolve(other: Path): Path = {
    if (other.isAbsolute) other
    else new NioGcsPath(chunks ++ other.asGcsPath.chunks, isAbsolute)
  }

  override def resolve(other: String): Path = {
    val otherPath = gcsFileSystem.getPath(other)
    if (otherPath.isAbsolute) otherPath
    else new NioGcsPath(chunks ++ otherPath.asGcsPath.chunks, isAbsolute)
  }

  override def toRealPath(options: LinkOption*): Path = this

  override def startsWith(other: Path): Boolean = {
    other match {
      case rel: NioGcsPath if !isAbsolute && rel.isAbsolute => false
      case _: NioGcsPath => chunks.startsWith(other.asGcsPath.chunks)
      case _ => false
    }
  }

  override def startsWith(other: String): Boolean = {
    Try(gcsFileSystem.getPath(other)) map {
      case rel: NioGcsPath if !isAbsolute && rel.isAbsolute => false
      case path@(_: NioGcsPath) => chunks.startsWith(path.asGcsPath.chunks)
      case _ => false
    } getOrElse false
  }

  override def toString: String = {
    if (absolute) s"$Protocol$fullPath"
    else fullPath
  }

  override def isAbsolute: Boolean = absolute
}
