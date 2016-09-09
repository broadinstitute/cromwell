package cromwell.util

import java.nio.file.Path

import better.files._
import wdl4s.values._

trait TestFileUtil {
  def createCannedFile(prefix: String, contents: String, dir: Option[Path] = None): Path = {
    val suffix = ".out"
    File.newTemporaryFile(prefix, suffix, dir.map(File.apply)).write(contents).path
  }

  def createFile(name: String, dir: Path, contents: String): Path = {
    File(dir).createDirectories()./(name).write(contents).path
  }
}

trait HashUtil extends TestFileUtil {
  // Files
  val file1 = WdlFile(createCannedFile("refFile", "some content").toString)
  val sameAsfile1 = WdlFile(createCannedFile("sameContent", "some content").toString)
  val anotherFile = WdlFile(createCannedFile("differentContent", "different content").toString)

  // Strings
  val string1 = WdlString("some text")
  val sameAsString1 = WdlString("some text")
  val anotherString = WdlString("different text")
}
