package cromwell.core.path

import java.io.IOException

import scala.util.{Failure, Try}

object PathCopier {

  /*
   * Remove common prefix found in path1 from path2 as long the directories match.
   */
  private[path] def truncateCommonRoot(path1: Path, path2: Path): String = {
    val string1 = path1.toAbsolutePath.pathWithoutScheme
    val string2 = path2.toAbsolutePath.pathWithoutScheme

    val regexIncludingSlashes = "(?<=/)|(?=/)" // http://stackoverflow.com/q/2206378
    val tokens1 = string1.split(regexIncludingSlashes)
    val tokens2 = string2.split(regexIncludingSlashes)

    val matchingTokens: Array[(String, String)] = tokens1.zip(tokens2).takeWhile(Function.tupled(_ == _))
    val matchingPrefix = matchingTokens.map({ case (str, _) => str }).mkString

    string2.stripPrefix(matchingPrefix).replaceAll("^/+", "")
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
    Option(destinationFilePath.parent).foreach(_.createPermissionedDirectories())
    Try {
      sourceFilePath.copyTo(destinationFilePath, overwrite = true)

      ()
    } recoverWith {
      case ex => Failure(new IOException(s"Failed to copy $sourceFilePath to $destinationFilePath", ex))
    }
  }
}
