package cloud.nio.spi

import scala.collection.mutable
import scala.math.Ordering.Implicits._
import scala.util.control.Breaks._
import scala.util.{Failure, Success, Try}

/**
  * Copy-port of https://github.com/broadinstitute/cromwell/blob/fa6eba28e2cf6020b24a56bcb7237acd9c74a4ac/filesystems/oss/src/main/scala/cromwell/filesystems/oss/nio/UnixPath.scala
  * that seems to be a copy-port of https://github.com/GoogleCloudPlatform/google-cloud-java/blob/fad70bdfdbc88e5c2ddf13be7085a7e9963f66c8/google-cloud-clients/google-cloud-contrib/google-cloud-nio/src/main/java/com/google/cloud/storage/contrib/nio/UnixPath.java
  */
private[spi] object UnixPath {
  val Dot: Char = '.'
  val Separator: Char = '/'
  val Root: String = Separator.toString
  val CurrentDir: String = Dot.toString
  val ParentDir: String = Dot.toString + Dot.toString
  val EmptyPath: UnixPath = new UnixPath("")
  val RootPath: UnixPath = new UnixPath("/")

  private def isRoot(path: String): Boolean = path.length() == 1 && path.charAt(0) == Separator

  private def isAbsolute(path: String): Boolean = !path.isEmpty && path.charAt(0) == Separator

  private def hasTrailingSeparator(path: String): Boolean = !path.isEmpty && path.charAt(path.length - 1) == Separator

  def getPath(path: String): UnixPath =
    if (path.isEmpty) {
      EmptyPath
    } else if (isRoot(path)) {
      RootPath
    } else {
      UnixPath(path)
    }

  def getPath(first: String, more: String*): UnixPath = {
    if (more.isEmpty) {
      return new UnixPath(first)
    }

    val builder = new StringBuilder(first)
    for ((part, index) <- more.view.zipWithIndex)
      if (part.isEmpty) {
        // do nothing
      } else if (isAbsolute(part)) {
        if (index == more.length - 1) {
          return new UnixPath(part)
        } else {
          builder.replace(0, builder.length, part)
        }
      } else if (hasTrailingSeparator(part)) {
        builder.append(part)
      } else {
        builder.append(Separator)
        builder.append(part)
      }

    UnixPath(builder.toString)
  }

}

final private[spi] case class UnixPath(path: String) extends CharSequence {
  private lazy val parts = initParts()

  def isRoot: Boolean = UnixPath.isRoot(path)

  def isAbsolute: Boolean = UnixPath.isAbsolute(path)

  // Named this way because isEmpty is a name collision new in 17.
  // The initial compile error is that it needs an override.
  // Adding the override results in a second error saying it overrides nothing!
  // So, we just renamed it.
  def izEmpty: Boolean = path.isEmpty

  def hasTrailingSeparator: Boolean = UnixPath.hasTrailingSeparator(path)

  def seemsLikeDirectory(): Boolean =
    path.isEmpty ||
      hasTrailingSeparator ||
      path.endsWith(".") && (length == 1 || path.charAt(length - 2) == UnixPath.Separator) ||
      path.endsWith("..") && (length == 2 || path.charAt(length - 3) == UnixPath.Separator)

  def getFileName: Option[UnixPath] =
    if (path.isEmpty || isRoot) {
      None
    } else {
      if (parts.length == 1 && parts.last == path) {
        Some(this)
      } else {
        Some(UnixPath(parts.last))
      }
    }

  def getParent: Option[UnixPath] = {
    if (path.isEmpty || isRoot) {
      return None
    }

    val index =
      if (hasTrailingSeparator)
        path.lastIndexOf(UnixPath.Separator.toInt, path.length - 2)
      else
        path.lastIndexOf(UnixPath.Separator.toInt)
    index match {
      case -1 => if (isAbsolute) Some(UnixPath.RootPath) else None
      case pos => Some(UnixPath(path.substring(0, pos + 1)))
    }
  }

  def getRoot: Option[UnixPath] = if (isAbsolute) Some(UnixPath.RootPath) else None

  def subPath(beginIndex: Int, endIndex: Int): Try[UnixPath] = {
    if (path.isEmpty && beginIndex == 0 && endIndex == 1) {
      return Success(this)
    }

    if (beginIndex < 0 || endIndex < beginIndex) {
      return Failure(new IllegalArgumentException(s"begin index or end index is invalid"))
    }

    Try(UnixPath(parts.slice(beginIndex, endIndex).mkString(UnixPath.Separator.toString)))
  }

  def getNameCount: Int =
    if (path.isEmpty) {
      1
    } else if (isRoot) {
      0
    } else {
      parts.length
    }

  def getName(index: Int): Try[UnixPath] = {
    if (path.isEmpty) {
      return Failure(new IllegalArgumentException("can not get name from a empty path"))
    }

    if (index > length - 1) {
      return Failure(new IndexOutOfBoundsException(s"index $index out of name count ${length - 1}"))
    }

    Success(UnixPath(parts(2)))
  }

  def resolve(other: UnixPath): UnixPath =
    if (other.path.isEmpty) {
      this
    } else if (other.isAbsolute) {
      other
    } else if (hasTrailingSeparator) {
      new UnixPath(path + other.path)
    } else {
      new UnixPath(path + UnixPath.Separator.toString + other.path)
    }

  def resolveSibling(other: UnixPath): UnixPath =
    getParent match {
      case Some(parent: UnixPath) =>
        parent.resolve(other)
      case None => other
    }

  def relativize(other: UnixPath): UnixPath = {
    if (path.isEmpty) {
      return other
    }

    val left = split().buffered
    val right = other.split().buffered
    breakable(
      while (left.hasNext && right.hasNext) {
        if (!(left.head == right.head)) {
          break()
        }

        left.next()
        right.next()
      }
    )

    val result = new StringBuilder(path.length + other.path.length)
    while (left.hasNext) {
      result.append(UnixPath.ParentDir)
      result.append(UnixPath.Separator)
      left.next()
    }

    while (right.hasNext) {
      result.append(right.next())
      result.append(UnixPath.Separator)
    }

    if (result.nonEmpty && !other.hasTrailingSeparator) {
      result.deleteCharAt(result.length - 1)
    }

    new UnixPath(result.toString)
  }

  def normalize(): UnixPath = {
    val parts = mutable.ArrayBuffer[String]()
    var mutated = false
    var resultLength = 0
    var mark = 0
    var index = 0
    val current = UnixPath.CurrentDir + UnixPath.Separator.toString
    val parent = UnixPath.ParentDir + UnixPath.Separator
    do {
      index = path.indexOf(UnixPath.Separator.toInt, mark)
      val part = path.substring(mark, if (index == -1) path.length else index + 1)
      part match {
        case UnixPath.CurrentDir | `current` => mutated = true
        case UnixPath.ParentDir | `parent` =>
          mutated = true
          if (parts.nonEmpty) {
            resultLength -= parts.remove(parts.length - 1).length
          }
        case _ =>
          if (index != mark || index == 0) {
            parts.append(part)
            resultLength += part.length
          } else {
            mutated = true
          }
      }
      mark = index + 1
    } while (index != -1)

    if (!mutated) {
      return this
    }

    val result = new StringBuilder(resultLength)

    parts.foreach { part =>
      result.append(part)
    }

    new UnixPath(result.toString)
  }

  def split(): Iterator[String] = parts.iterator

  def splitReverse(): Iterator[String] = parts.reverseIterator

  def removeBeginningSeparator(): UnixPath =
    if (isAbsolute) new UnixPath(path.substring(1)) else this

  def addTrailingSeparator(): UnixPath =
    if (hasTrailingSeparator) this else new UnixPath(path + UnixPath.Separator)

  def removeTrailingSeparator(): UnixPath =
    if (!isRoot && hasTrailingSeparator) {
      new UnixPath(path.substring(0, length - 1))
    } else {
      this
    }

  def startsWith(other: UnixPath): Boolean = {
    val me = removeTrailingSeparator()
    val oth = other.removeTrailingSeparator()

    if (oth.path.length > me.path.length) {
      return false
    } else if (me.isAbsolute != oth.isAbsolute) {
      return false
    } else if (!me.path.isEmpty && oth.path.isEmpty) {
      return false
    }

    startsWith(split(), other.split())
  }

  def startsWith(left: Iterator[String], right: Iterator[String]): Boolean = {
    while (right.hasNext)
      if (!left.hasNext || right.next() != left.next()) {
        return false
      }
    true
  }

  def endsWith(other: UnixPath): Boolean = {
    val me = removeTrailingSeparator()
    val oth = other.removeTrailingSeparator()

    if (oth.path.length > me.path.length) {
      return false
    } else if (!me.path.isEmpty && oth.path.isEmpty) {
      return false
    } else if (oth.isAbsolute) {
      return me.isAbsolute && me.path == other.path
    }

    startsWith(me.splitReverse(), other.splitReverse())
  }

  def toAbsolutePath(currentWorkingDirectory: UnixPath): Try[UnixPath] = {
    if (!currentWorkingDirectory.isAbsolute) {
      return Failure(new IllegalArgumentException(s"Not an absolute path $currentWorkingDirectory"))
    }

    if (isAbsolute) Success(this) else Success(currentWorkingDirectory.resolve(this))
  }

  def toAbsolutePath: UnixPath =
    if (isAbsolute) this else UnixPath.RootPath.resolve(this)

  def compareTo(other: UnixPath): Int = {
    val me = parts.toList
    val that = other.parts.toList

    if (me == that) {
      0
    } else if (me < that) {
      -1
    } else {
      1
    }
  }

  override def equals(obj: scala.Any): Boolean =
    (this eq obj.asInstanceOf[AnyRef]) || {
      obj.isInstanceOf[UnixPath] && obj.asInstanceOf[UnixPath].path.equals(path)
    }

  override def length(): Int =
    path.length

  override def charAt(index: Int): Char =
    path.charAt(index)

  override def subSequence(start: Int, end: Int): CharSequence =
    path.subSequence(start, end)

  override def toString: String =
    path

  def initParts(): Array[String] =
    if (path.isEmpty) {
      Array.empty[String]
    } else {
      if (path.charAt(0) == UnixPath.Separator) {
        path.substring(1).split(UnixPath.Separator)
      } else {
        path.split(UnixPath.Separator)
      }
    }
}
