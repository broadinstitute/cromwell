package cromwell.core.path

import java.io.IOException
import java.nio.file.Path

import better.files._

import scala.util.{Failure, Try}

object PathCopier {
  def getDestinationFilePath(sourceContextPath: Path, sourceFilePath: Path, destinationDirPath: Path): Path = {
    val relativeFileString = sourceContextPath.toAbsolutePath.relativize(sourceFilePath.toAbsolutePath).toString
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
