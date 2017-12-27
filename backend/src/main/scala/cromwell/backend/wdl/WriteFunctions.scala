package cromwell.backend.wdl

import better.files.File.OpenOptions
import cats.instances.future._
import cats.syntax.functor._
import cromwell.core.io.AsyncIoFunctions
import cromwell.core.path.Path
import wom.expression.IoFunctionSet
import wom.values.WomSingleFile

import scala.concurrent.Future

trait WriteFunctions extends IoFunctionSet with AsyncIoFunctions {

  /**
    * Directory that will be used to write files.
    */
  def writeDirectory: Path

  private lazy val _writeDirectory = writeDirectory.createPermissionedDirectories()

  override def writeFile(path: String, content: String): Future[WomSingleFile] = {
    val file = _writeDirectory / path
    asyncIo.existsAsync(file) flatMap {
      case false => asyncIo.writeAsync(file, content, OpenOptions.default) as { WomSingleFile(file.pathAsString) }
      case true => Future.successful(WomSingleFile(file.pathAsString))
    }
  }
}
