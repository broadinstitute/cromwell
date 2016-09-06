package wdl4s

import java.nio.file.{Files, Path}

import better.files._

trait TestFileUtil {
  def createCannedFile(prefix: String, contents: String, dir: Option[Path] = None): Path = {
    val suffix = ".out"
    val file = dir match {
      case Some(path) => Files.createTempFile(path, prefix, suffix)
      case None => Files.createTempFile(prefix, suffix)
    }
    File(file).write(contents).path
  }

  def createFile(name: String, dir: Path, contents: String) = {
    File(dir).createDirectories()./(name).write(contents).path
  }
}
