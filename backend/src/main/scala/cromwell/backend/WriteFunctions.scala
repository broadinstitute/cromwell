package cromwell.backend

import better.files.File.OpenOptions
import cats.implicits._
import cromwell.core.io.AsyncIoFunctions
import cromwell.core.path.{Path, PathFactory}
import wom.expression.IoFunctionSet
import wom.values.WomSingleFile

import scala.concurrent.Future

trait WriteFunctions extends PathFactory with IoFunctionSet with AsyncIoFunctions {

  /**
    * True if the writeDirectory will be mounted onto a docker container
    */
  def isDocker: Boolean

  /**
    * Directory that will be used to write files.
    */
  def writeDirectory: Path

  private lazy val _writeDirectory =
    if (isDocker) writeDirectory.createPermissionedDirectories() else writeDirectory.createDirectories()

  protected def writeAsync(file: Path, content: String) = asyncIo.writeAsync(file, content, OpenOptions.default)

  override def writeFile(path: String, content: String): Future[WomSingleFile] = {
    val file = _writeDirectory / path
    asyncIo.existsAsync(file) flatMap {
      case false => writeAsync(file, content) as WomSingleFile(file.pathAsString)
      case true => Future.successful(WomSingleFile(file.pathAsString))
    }
  }
}
