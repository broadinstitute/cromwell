package cromwell.backend.wdl

import better.files.File.OpenOptions
import cats.instances.future._
import cats.syntax.functor._
import common.util.StringUtil._
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

  private lazy val _writeDirectory = if (isDocker) writeDirectory.createPermissionedDirectories() else writeDirectory.createDirectories()

  override def writeFile(path: String, content: String): Future[WomSingleFile] = {
    val file = _writeDirectory / path
    asyncIo.existsAsync(file) flatMap {
      case false => asyncIo.writeAsync(file, content, OpenOptions.default) as { WomSingleFile(file.pathAsString) }
      case true => Future.successful(WomSingleFile(file.pathAsString))
    }
  }

  private val relativeToLocal = System.getProperty("user.dir").ensureSlashed

  def relativeToAbsolutePath(pathFrom: String): String = if (buildPath(pathFrom).isAbsolute) pathFrom else relativeToLocal + pathFrom

  override def copyFile(pathFrom: String, targetName: String): Future[WomSingleFile] = {
    val source = buildPath(relativeToAbsolutePath(pathFrom))
    val destination = _writeDirectory / targetName

    asyncIo.copyAsync(source, destination).as(WomSingleFile(destination.pathAsString)) recoverWith {
      case e => Future.failed(new Exception(s"Could not copy ${source.toAbsolutePath} to ${destination.toAbsolutePath}", e))
    }
  }
}
