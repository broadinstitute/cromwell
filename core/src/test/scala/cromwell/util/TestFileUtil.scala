package cromwell.util

import java.nio.file.attribute.PosixFilePermission

import cromwell.core.path.{DefaultPathBuilder, Path}
import wdl4s.wdl.values._

trait TestFileUtil {
  def createCannedFile(prefix: String, contents: String, dir: Option[Path] = None): Path = {
    val suffix = ".out"
    val tempFile = DefaultPathBuilder.createTempFile(prefix, suffix, dir)
    tempFile.createPermissionedDirectories()
    tempFile.addPermission(PosixFilePermission.OTHERS_READ)
    tempFile.addPermission(PosixFilePermission.OTHERS_WRITE)
    tempFile.write(contents)
  }

  def createFile(name: String, dir: Path, contents: String): Path = {
    dir.createPermissionedDirectories()./(name).write(contents)
  }
}

trait HashUtil extends TestFileUtil {
  // Files
  val file1 = WdlFile(createCannedFile("refFile", "some content").pathAsString)
  val sameAsfile1 = WdlFile(createCannedFile("sameContent", "some content").pathAsString)
  val anotherFile = WdlFile(createCannedFile("differentContent", "different content").pathAsString)

  // Strings
  val string1 = WdlString("some text")
  val sameAsString1 = WdlString("some text")
  val anotherString = WdlString("different text")
}
