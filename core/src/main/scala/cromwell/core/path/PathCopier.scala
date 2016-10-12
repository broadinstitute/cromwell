package cromwell.core.path

import java.io.IOException
import java.nio.file.Path

import better.files._

import scala.util.{Failure, Try}

object PathCopier {

  /*
   * Remove p1 from p2 as long as they match.
   */
  private def truncateCommonRoot(p1: Path, p2: Path): String = {
    def names(p: Path) = 0 until p.getNameCount map p.getName

    val names1 = names(p1)

    val truncated = names(p2).zipWithIndex.dropWhile {
      case (n1, n2) => n2 < names1.size && n1.equals(names1(n2))
    } map { _._1 }

    truncated match {
      case empty if empty.isEmpty => ""
      case truncs => truncs.reduceLeft(_.resolve(_)).toString
    }
  }

  def getDestinationFilePath(sourceContextPath: Path, sourceFilePath: Path, destinationDirPath: Path): Path = {
    val relativeFileString = truncateCommonRoot(sourceContextPath.toAbsolutePath, sourceFilePath.toAbsolutePath)
    destinationDirPath.resolve(relativeFileString)
  }

  /**
    * Copies from a relative source to destination dir. NOTE: Copies are not atomic, and may create a partial copy.
    */
  def copy(sourceContextPath: Path, sourceFilePath: Path, destinationDirPath: Path): Try[Unit] = {
    val destinationFilePath = getDestinationFilePath(sourceContextPath, sourceFilePath, destinationDirPath)
    copy(sourceFilePath, destinationFilePath)
  }

  /**
    * Copies from source to destination. NOTE: Copies are not atomic, and may create a partial copy.
    */
  def copy(sourceFilePath: Path, destinationFilePath: Path): Try[Unit] = {
    Option(File(destinationFilePath).parent).foreach(_.createDirectories())
    Try {
      File(sourceFilePath).copyTo(destinationFilePath, overwrite = true)

      ()
    } recoverWith {
      case ex => Failure(new IOException(s"Failed to copy ${sourceFilePath.toUri} to ${destinationFilePath.toUri}", ex))
    }
  }
}
