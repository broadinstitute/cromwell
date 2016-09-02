package cromwell.core

import java.io.Writer
import java.nio.file.{FileSystem, Path}

import scala.collection.immutable.Queue
import scala.util.{Success, Try}

trait PathFactory {
  def findFileSystem(rawString: String, fss: List[FileSystem], mapping: PartialFunction[FileSystem, Try[Path]]) = {
    fss.toStream collect mapping collectFirst { case Success(p) => p } getOrElse {
      throw new IllegalArgumentException(s"Could not find suitable filesystem to parse $rawString")
    }
  }

  def buildPath(rawString: String, fileSystems: List[FileSystem]): Path = {
    findFileSystem(rawString, fileSystems, {
      case fs: FileSystem => Try(fs.getPath(rawString))
    })
  }
}

object PathFactory {
  def swapExt(filePath: String, oldExt: String, newExt: String): String = {
    filePath.stripSuffix(oldExt) + newExt
  }

  implicit class EnhancedPath(val path: Path) extends AnyVal {
    def swapExt(oldExt: String, newExt: String): Path = {
      path.getFileSystem.getPath(s"${path.toString.stripSuffix(oldExt)}$newExt")
    }

    def untailed = UntailedWriter(path)

    def tailed(tailedSize: Int) = TailedWriter(path, tailedSize)
  }

  implicit class FlushingAndClosingWriter(writer: Writer) {
    /** Convenience method to flush and close in one shot. */
    def flushAndClose() = {
      writer.flush()
      writer.close()
    }
  }
}

/**
  * Used with a `ProcessLogger`, writes lines with a newline.
  */
trait PathWriter {
  import better.files._

  val path: Path
  lazy val writer: Writer = path.newBufferedWriter

  /**
    * Passed to `ProcessLogger` to add a new line.
    *
    * @param string Line to add to the logs.
    */
  def writeWithNewline(string: String) {
    writer.write(string)
    writer.write("\n")
  }
}

/**
  * Used with a `ProcessLogger`, writes lines with a newline.
  *
  * @param path Path to the log file.
  */
case class UntailedWriter(path: Path) extends PathWriter

/**
  * Used with a `ProcessLogger`, queues up the `tailedSize` number of lines.
  *
  * @param path Path to the log file.
  * @param tailedSize Maximum number of lines to save in the internal FIFO queue.
  */
case class TailedWriter(path: Path, tailedSize: Int) extends PathWriter {
  private var isTailed = false
  private var tailedLines: Queue[String] = Queue.empty

  /**
    * Passed to `ProcessLogger` to add a new line, and adds the line to the tailed queue.
    *
    * @param string Line to add to the logs.
    */
  override def writeWithNewline(string: String) {
    tailedLines :+= string
    while (tailedLines.size > tailedSize) {
      tailedLines = tailedLines.takeRight(tailedSize)
      isTailed = true
    }
    super.writeWithNewline(string)
  }

  /**
    * Returns a descriptive tail of the `path` and the last `tailedLines` written.
    *
    * @return a descriptive tail of the `path` and the last `tailedLines` written.
    */
  def tailString: String = {
    if (tailedLines.isEmpty) {
      s"Contents of $path were empty."
    } else if (isTailed) {
      s"Last ${tailedLines.size} of $path:\n${tailedLines.mkString("\n")}"
    } else {
      s"Contents of $path:\n${tailedLines.mkString("\n")}"
    }
  }
}

