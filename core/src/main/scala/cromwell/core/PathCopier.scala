package cromwell.core

import java.io.IOException
import java.nio.file.Path

import better.files._

object PathCopier {
  def getDestinationFilePath(sourceContextPath: Path, sourceFilePath: Path, destinationDirPath: Path): Path = {
    val relativeFileString = sourceContextPath.toAbsolutePath.relativize(sourceFilePath.toAbsolutePath).toString
    destinationDirPath.resolve(relativeFileString)
  }

  /**
    * Copies from a relative source to destination dir. NOTE: Copies are not atomic, and may create a partial copy.
    */
  def copy(sourceContextPath: Path, sourceFilePath: Path, destinationDirPath: Path): Unit = {
    val destinationFilePath = getDestinationFilePath(sourceContextPath, sourceFilePath, destinationDirPath)
    copy(sourceFilePath, destinationFilePath)
  }

  /**
    * Copies from source to destination. NOTE: Copies are not atomic, and may create a partial copy.
    */
  def copy(sourceFilePath: Path, destinationFilePath: Path): Unit = {
    Option(File(destinationFilePath).parent).foreach(_.createDirectories())
    try {
      File(sourceFilePath).copyTo(destinationFilePath, overwrite = true)
      ()
    } catch {
      case ex: Exception =>
        throw new IOException(s"Failed to copy ${sourceFilePath.toUri.toString} to ${destinationFilePath.toUri.toString}", ex)
    }
  }
}
