package cromwell.util

import java.io.{File, FileInputStream, Writer}
import java.nio.file.Path

import better.files._
import org.apache.commons.codec.digest.DigestUtils
import wdl4s.values.Hashable

import scala.collection.immutable.Queue
import scala.util.{Failure, Success, Try}

object FileUtil {
  def swapExt(filePath: String, oldExt: String, newExt: String): String = {
    filePath.stripSuffix(oldExt) + newExt
  }

  def parseTsv(tsv: String): Try[Array[Array[String]]] = {
    val table = tsv.split("\n").map(_.split("\t"))
    table.map(_.size).toSet match {
      case s if s.size > 1 => Failure(new UnsupportedOperationException("TSV is not uniform"))
      case _ => Success(table)
    }
  }

  implicit class FlushingAndClosingWriter(writer: Writer) {
    /** Convenience method to flush and close in one shot. */
    def flushAndClose() = {
      writer.flush()
      writer.close()
    }
  }

  implicit class EnhancedFile(val file: File) extends AnyVal with Hashable {
    /** Read an entire file into a string, closing the underlying stream. */
    def slurp: String = {
      // TODO: deprecate slurp, and java.io.File in general?
      file.toPath.contentAsString
    }

    def md5Sum: String = {
      val fis = new FileInputStream(file)
      try {
        DigestUtils.md5Hex(fis)
      } finally fis.close()
    }
  }

  implicit class EnhancedPath(val path: Path) extends AnyVal {
    def swapExt(oldExt: String, newExt: String): Path = {
      path.getFileSystem.getPath(FileUtil.swapExt(path.toString, oldExt, newExt))
    }

    def untailed = new UntailedWriter(path)

    def tailed(tailedSize: Int) = TailedWriter(path, tailedSize)
  }
}

/**
 * Used with a `ProcessLogger`, writes lines with a newline.
 */
trait PathWriter {
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
  var isTailed = false
  var tailedLines: Queue[String] = Queue.empty

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
