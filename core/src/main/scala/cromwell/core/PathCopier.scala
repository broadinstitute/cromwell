package cromwell.core

import java.nio.file.Path

import better.files._

object PathCopier {
  def copy(sourceContextPath: Path, sourceFilePath: Path, destinationDirPath: Path, overwrite: Boolean): Unit = {
    val relativeFileString: String = sourceContextPath.toAbsolutePath.relativize(sourceFilePath.toAbsolutePath).toString
    val destinationFilePath: Path = destinationDirPath.resolve(relativeFileString)
    copy(sourceFilePath, destinationFilePath, overwrite)
  }

  def copy(sourceFilePath: Path, destinationFilePath: Path, overwrite: Boolean = false): Unit = {
    Option(File(destinationFilePath).parent).foreach(_.createDirectories())
    File(sourceFilePath).copyTo(destinationFilePath, overwrite)
  }
}
