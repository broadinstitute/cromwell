package centaur.test.workflow

import java.nio.file.Path
import java.nio.file.Files

import better.files.File
import cats.data.Validated._
import centaur.test._
import com.typesafe.config.Config
import configs.Result
import configs.syntax._
import spray.json.JsArray

case class WorkflowData(wdl: String, inputs: Option[String], options: Option[String], zippedImports: Option[File])

object WorkflowData {
  def fromConfig(conf: Config, basePath: Path): ErrorOr[WorkflowData] = {
    conf.get[Path]("wdl") match {
      case Result.Success(wdl) => Valid(WorkflowData(basePath.resolve(wdl), conf, basePath))
      case Result.Failure(_) => invalidNel("No wdl path provided")
    }
  }

  def apply(wdl: Path, conf: Config, basePath: Path): WorkflowData = {
    def getOptionalPath(name: String) = conf.get[Option[Path]](name) valueOrElse None map basePath.resolve

    def getImports = conf.get[List[Path]]("imports") match {
      case Result.Success(paths) => zipImports(paths map basePath.resolve)
      case Result.Failure(_) => None
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

    // TODO: The slurps can throw - not a high priority but see #36
    WorkflowData(wdl.slurp, getOptionalPath("inputs") map { _.slurp }, getOptionalPath("options") map { _.slurp }, getImports)
  }

  implicit class EnhancedPath(val path: Path) extends AnyVal {
    /** Read an entire file into a string, closing the underlying stream. */
    def slurp: String = {
      val source = io.Source.fromFile(path.toFile, "UTF-8")
      try source.mkString finally source.close()
    }
  }
}
