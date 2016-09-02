package cromwell.util

import java.nio.file.{Files, Path}

import better.files._
import wdl4s.values._

trait TestFileUtil {
  /* TODO: Copy pasted from wdl4s's test artifact. */
  def createCannedFile(prefix: String, contents: String, dir: Option[Path] = None): Path = {
    val suffix = ".out"
    val file = dir match {
      case Some(path) => Files.createTempFile(path, prefix, suffix)
      case None => Files.createTempFile(prefix, suffix)
    }
    file.write(contents).path
  }

  /* TODO: Copy pasted from wdl4s's test artifact. */
  def createFile(name: String, dir: Path, contents: String): Path = {
    dir.createDirectories()./(name).write(contents).path
  }
}

trait HashUtil extends TestFileUtil {
  // Files
  val file1 = WdlFile(createCannedFile("refFile", "some content").fullPath)
  val sameAsfile1 = WdlFile(createCannedFile("sameContent", "some content").fullPath)
  val anotherFile = WdlFile(createCannedFile("differentContent", "different content").fullPath)

  // Strings
  val string1 = WdlString("some text")
  val sameAsString1 = WdlString("some text")
  val anotherString = WdlString("different text")
}
