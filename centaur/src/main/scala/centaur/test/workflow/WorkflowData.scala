package centaur.test.workflow

import java.nio.file.Path

import better.files.File
import cats.data.Validated._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import configs.Result.{Failure, Success}
import configs.syntax._
import cromwell.api.model.Label
import net.ceedubs.ficus.Ficus._
import spray.json._


case class WorkflowData(workflowContent: Option[String],
                        workflowUrl: Option[String],
                        workflowRoot: Option[String],
                        workflowType: Option[String],
                        workflowTypeVersion: Option[String],
                        inputs: Option[String],
                        options: Option[String],
                        labels: List[Label],
                        zippedImports: Option[File])

object WorkflowData {
  def fromConfig(filesConfig: Config, fullConfig: Config, basePath: Path): ErrorOr[WorkflowData] = {
    val workflowUrl = filesConfig.as[Option[String]]("workflowUrl")
    val workflowSourcePath = filesConfig.as[Option[String]]("workflow")

    (workflowSourcePath, workflowUrl) match {
      case (Some(workflowPath), None) => Valid(WorkflowData(
        workflowPath = Option(basePath.resolve(workflowPath)),
        workflowUrl = None,
        filesConfig = filesConfig,
        fullConfig = fullConfig,
        basePath = basePath))
      case (None, Some(_)) => Valid(WorkflowData(
        workflowPath = None,
        workflowUrl = workflowUrl,
        filesConfig = filesConfig,
        fullConfig = fullConfig,
        basePath = basePath))
      case (Some(_), Some(_)) => invalidNel(s"Both 'workflow' path or 'workflowUrl' can't be provided.")
      case (None, None) => invalidNel(s"No 'workflow' path or 'workflowUrl' provided.")
    }
  }

  def apply(workflowPath: Option[Path], workflowUrl: Option[String], filesConfig: Config, fullConfig: Config, basePath: Path): WorkflowData = {
    def getOptionalPath(name: String) = filesConfig.get[Option[Path]](name) valueOrElse None map basePath.resolve

    def getImports = filesConfig.get[List[Path]]("imports") match {
      case Success(paths) => zipImports(paths map basePath.resolve)
      case Failure(_) => None
    }

    def zipImports(imports: List[Path]): Option[File] = {
      val zippedDir = imports match {
        case Nil => None
        case _ =>
          val importsDirName = workflowPath.get.getFileName.toString.replaceAll("\\.[^.]*$", "") //TODO: Saloni-What about this?
          val importsDir = File.newTemporaryDirectory(importsDirName + "_imports")
          imports foreach { p =>
            val srcFile = File(p.toAbsolutePath.toString)
            val destFile = importsDir / srcFile.name
            srcFile.copyTo(destFile, overwrite = true) }

          Option(importsDir.zip())
      }
      zippedDir
    }

    def getLabels: List[Label] = {
      import cromwell.api.model.LabelsJsonFormatter._
      getOptionalPath("labels") map { _.slurp.parseJson.convertTo[List[Label]] } getOrElse List.empty
    }

    val workflowSource = if (workflowPath.isDefined) Option(workflowPath.get.slurp) else None
    val workflowType = fullConfig.get[Option[String]]("workflowType").value
    val workflowTypeVersion = fullConfig.get[Option[String]]("workflowTypeVersion").value
    val workflowRoot = fullConfig.get[Option[String]]("workflowRoot").value

    // TODO: The slurps can throw - not a high priority but see #36
    WorkflowData(
      workflowContent = workflowSource,
      workflowUrl = workflowUrl,
      workflowRoot = workflowRoot,
      workflowType = workflowType,
      workflowTypeVersion = workflowTypeVersion,
      inputs = getOptionalPath("inputs") map { _.slurp },
      options = getOptionalPath("options") map { _.slurp },
      labels = getLabels,
      zippedImports = getImports
    )
  }

  implicit class EnhancedPath(val path: Path) extends AnyVal {
    /** Read an entire file into a string, closing the underlying stream. */
    def slurp: String = {
      val source = scala.io.Source.fromFile(path.toFile, "UTF-8")
      try source.mkString finally source.close()
    }
  }
}
