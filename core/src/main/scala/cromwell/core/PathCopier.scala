package cromwell.core

import java.nio.file.Path

import better.files._

object PathCopier {
  def copy(sourceContextPath: Path, sourceFilePath: Path, destinationDirPath: Path): Unit = {
    val relativeFileString: String = sourceContextPath.toAbsolutePath.relativize(sourceFilePath.toAbsolutePath).toString
    val destinationFilePath: Path = destinationDirPath.resolve(relativeFileString)
    copy(sourceFilePath, destinationFilePath)
  }

  def copy(sourceFilePath: Path, destinationFilePath: Path): Unit = {
    Option(destinationFilePath.parent).foreach(_.createDirectories())
    sourceFilePath.copyTo(destinationFilePath)
  }
}
