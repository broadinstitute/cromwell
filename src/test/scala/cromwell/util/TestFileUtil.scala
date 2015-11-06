package cromwell.util

import java.io.{FileWriter, File}
import java.nio.file.{Files, Path}

import cromwell.binding.values._

trait TestFileUtil {
  private def write(file: File, contents: String) = {
    val writer = new FileWriter(file)
    writer.write(contents)
    writer.flush()
    writer.close()
    file
  }

  def createCannedFile(prefix: String, contents: String, dir: Option[Path] = None): File = {
    val suffix = ".out"
    val file = dir match {
      case Some(path) => Files.createTempFile(path, prefix, suffix)
      case None => Files.createTempFile(prefix, suffix)
    }
    write(file.toFile, contents)
  }

  def createFile(name: String, dir: Path, contents: String) = {
    dir.toFile.mkdirs()
    write(dir.resolve(name).toFile, contents)
  }
}

trait HashUtil extends TestFileUtil {
  // Files
  val file1 = WdlFile(createCannedFile("refFile", "some content").toPath.toAbsolutePath.toString)
  val sameAsfile1 = WdlFile(createCannedFile("sameContent", "some content").toPath.toAbsolutePath.toString)
  val anotherFile = WdlFile(createCannedFile("differentContent", "different content").toPath.toAbsolutePath.toString)

  // Strings
  val string1 = WdlString("some text")
  val sameAsString1 = WdlString("some text")
  val anotherString = WdlString("different text")

}