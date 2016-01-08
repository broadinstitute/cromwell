package wdl4s

import java.io.{FileWriter, File}
import java.nio.file.{Files, Path}

import wdl4s.values._

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

