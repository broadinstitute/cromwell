package centaur.test.workflow

import java.nio.file.Path

import better.files.File
import cats.data.Validated._
import centaur.test._
import com.typesafe.config.Config
import configs.Result.{Failure, Success}
import configs.syntax._
import cromwell.api.model.Label
import spray.json._


case class WorkflowData(wdl: String, inputs: Option[String], options: Option[String], labels: List[Label], zippedImports: Option[File])

object WorkflowData {
  def fromConfig(conf: Config, basePath: Path): ErrorOr[WorkflowData] = {
    conf.get[Path]("wdl") match {
      case Success(wdl) => Valid(WorkflowData(basePath.resolve(wdl), conf, basePath))
      case Failure(_) => invalidNel("No wdl path provided")
    }
  }

  def apply(wdl: Path, conf: Config, basePath: Path): WorkflowData = {
    def getOptionalPath(name: String) = conf.get[Option[Path]](name) valueOrElse None map basePath.resolve

    def getImports = conf.get[List[Path]]("imports") match {
      case Success(paths) => zipImports(paths map basePath.resolve)
      case Failure(_) => None
    }

    def zipImports(imports: List[Path]): Option[File] = {
      val zippedDir = imports match {
        case x :: _ => {
          val importsDirName = wdl.getFileName.toString.replaceAll("\\.[^.]*$", "")
          val importsDir = File.newTemporaryDirectory(importsDirName + "_imports")
          imports foreach { p =>
            val srcFile = File(p.toAbsolutePath.toString)
            val destFile = importsDir / srcFile.name
            srcFile.copyTo(destFile, true) }

          Option(importsDir.zip())
        }
        case Nil => None
      }
      zippedDir
    }

    def getLabels: List[Label] = {
      import cromwell.api.model.LabelsJsonFormatter._
      getOptionalPath("labels") map { _.slurp } map { _.parseJson.convertTo[List[Label]] } getOrElse List.empty
    }

    // TODO: The slurps can throw - not a high priority but see #36
    WorkflowData(wdl.slurp, getOptionalPath("inputs") map { _.slurp }, getOptionalPath("options") map { _.slurp }, getLabels, getImports)
  }

  implicit class EnhancedPath(val path: Path) extends AnyVal {
    /** Read an entire file into a string, closing the underlying stream. */
    def slurp: String = {
      val source = scala.io.Source.fromFile(path.toFile, "UTF-8")
      try source.mkString finally source.close()
    }
  }
}
