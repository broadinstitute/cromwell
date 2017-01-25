package cromwell.util

import java.nio.file.Path
import java.nio.file.attribute.PosixFilePermission

import better.files._
import wdl4s.values._
import cromwell.core.path.FileImplicits._

trait TestFileUtil {
  def createCannedFile(prefix: String, contents: String, dir: Option[Path] = None): Path = {
    val suffix = ".out"
    val tempFile = File.newTemporaryFile(prefix, suffix, dir.map(File.apply))
    tempFile.createPermissionedDirectories()
    tempFile.addPermission(PosixFilePermission.OTHERS_READ)
    tempFile.addPermission(PosixFilePermission.OTHERS_WRITE)
    tempFile.write(contents).path
  }

  def createFile(name: String, dir: Path, contents: String): Path = {
    File(dir).createPermissionedDirectories()./(name).write(contents).path
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
