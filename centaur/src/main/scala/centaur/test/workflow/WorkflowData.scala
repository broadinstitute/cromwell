package centaur.test.workflow

import better.files.File
import cats.data.Validated._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import configs.Result.{Failure, Success}
import configs.syntax._
import cromwell.api.model.Label
import net.ceedubs.ficus.Ficus._
import spray.json._

/**
  * Stores workflow data / references.
  *
  * For some test suites the inputs and options may dynamically rendered with secure variables. For the other test
  * suites that don't require these inputs / options, the respective Function1[Unit, String] will not be invoked to read
  * the non-existent content.
  */
case class WorkflowData(workflowContent: Option[String],
                        workflowUrl: Option[String],
                        workflowRoot: Option[String],
                        workflowType: Option[String],
                        workflowTypeVersion: Option[String],
                        inputs: Option[() => String],
                        options: Option[() => String],
                        labels: List[Label],
                        zippedImports: Option[File])

object WorkflowData {
  def fromConfig(filesConfig: Config, fullConfig: Config, basePath: File): ErrorOr[WorkflowData] = {
    val workflowUrl = filesConfig.as[Option[String]]("workflowUrl")
    val workflowSourcePath = filesConfig.as[Option[String]]("workflow")

    (workflowSourcePath, workflowUrl) match {
      case (Some(workflowPath), None) => Valid(WorkflowData(
        workflowPath = Option(basePath / workflowPath),
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

  def apply(workflowPath: Option[File], workflowUrl: Option[String], filesConfig: Config, fullConfig: Config, basePath: File): WorkflowData = {
    def getOptionalFile(name: String): Option[File] = {
      filesConfig.get[Option[String]](name) valueOrElse None map basePath./
    }

    def getImports = filesConfig.get[List[String]]("imports") match {
      case Success(paths) => zipImports(paths map basePath./)
      case Failure(_) => None
    }

    def getImportsDirName(workflowPath: Option[File], workflowUrl: Option[String]): String = {
      workflowPath match {
        case Some(file) => file.name.replaceAll("\\.[^.]*$", "")
        case None => // workflow url is defined
          val fileName = workflowUrl.get.split("/").last
          fileName.replaceAll("\\.[^.]*$", "")
      }
    }

    def zipImports(imports: List[File]): Option[File] = {
      imports match {
        case Nil => None
        case _ =>
          val importsDirName = getImportsDirName(workflowPath, workflowUrl)
          val importsDir = File.newTemporaryDirectory(importsDirName + "_imports")
          imports foreach { p =>
            val srcFile = File(p.pathAsString)
            val destFile = importsDir / srcFile.name
            srcFile.copyTo(destFile, overwrite = true)
          }

          Option(importsDir.zip())
      }
    }

    def getLabels: List[Label] = {
      import cromwell.api.model.LabelsJsonFormatter._
      getOptionalFile("labels") map { _.contentAsString.parseJson.convertTo[List[Label]] } getOrElse List.empty
    }

    val workflowSource = if (workflowPath.isDefined) Option(workflowPath.get.contentAsString) else None
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
      inputs = getOptionalFile("inputs") map (file => () => file.contentAsString),
      options = getOptionalFile("options") map (file => () => file.contentAsString),
      labels = getLabels,
      zippedImports = getImports
    )
  }
}
