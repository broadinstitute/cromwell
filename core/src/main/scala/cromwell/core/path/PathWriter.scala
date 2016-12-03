package cromwell.core.path

import java.io.Writer
import java.nio.file.Path

import scala.collection.immutable.Queue



/**
  * Used with a `ProcessLogger`, writes lines with a newline.
  */
trait PathWriter {
  import better.files._

  val path: Path
  lazy val writer: Writer = File(path).newBufferedWriter

  /**
    * Passed to `ProcessLogger` to add a new line.
    *
    * @param string Line to add to the logs.
    */
  def writeWithNewline(string: String): Unit = {
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
  override def writeWithNewline(string: String): Unit = {
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

