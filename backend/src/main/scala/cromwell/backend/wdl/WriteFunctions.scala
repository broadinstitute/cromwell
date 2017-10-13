package cromwell.backend.wdl

import cromwell.core.path.Path
import wom.expression.IoFunctionSet
import wom.values.WdlFile

import scala.concurrent.Future
import scala.util.Try

trait WriteFunctions extends IoFunctionSet {

  /**
    * Directory that will be used to write files.
    */
  def writeDirectory: Path

  private lazy val _writeDirectory = writeDirectory.createPermissionedDirectories()

  override def writeFile(path: String, content: String): Future[WdlFile] = {
    val file = _writeDirectory / path
    Future.fromTry(
      Try(if (file.notExists) file.write(content)) map { _ => WdlFile(file.pathAsString) }
    )
  }
}
