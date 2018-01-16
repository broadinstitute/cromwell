package cromwell.backend.wdl

import java.io.File

import better.files.File.OpenOptions
import cats.instances.future._
import cats.syntax.functor._
import com.typesafe.config.ConfigFactory
import cromwell.core.io.AsyncIoFunctions
import cromwell.core.path.{Path, PathFactory}
import wom.expression.IoFunctionSet
import wom.values.WomSingleFile

import net.ceedubs.ficus.Ficus._
import scala.concurrent.Future


trait WriteFunctions extends PathFactory with IoFunctionSet with AsyncIoFunctions {

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

  private val relativeToLocal = System.getProperty("user.dir")
  private val relativeToPapi = ConfigFactory.load().as[Option[String]]("papi.default-input-gcs-prefix").orElse(ConfigFactory.load().as[Option[String]]("backend.providers.JES.config.root")).getOrElse("gs://")

  private def withSlash(str: String) = if (str.endsWith("/")) str else str + "/"

  def relativeToAbsolutePath(pathFrom: String, papi: Boolean): String = if (new File(pathFrom).isAbsolute) pathFrom else {
    val result = withSlash(if (papi) relativeToPapi else relativeToLocal) + pathFrom
    println(s"RelativeToAbsolute of $pathFrom was $result")
    result
  }

  override def copyFile(pathFrom: String, targetName: String): Future[WomSingleFile] = {
    val papi = buildPath(".").toString.startsWith("gs")

    val source = buildPath(relativeToAbsolutePath(pathFrom, papi))
    val destination = _writeDirectory / targetName

    asyncIo.copyAsync(source, destination).as(WomSingleFile(destination.pathAsString))
  }
}
