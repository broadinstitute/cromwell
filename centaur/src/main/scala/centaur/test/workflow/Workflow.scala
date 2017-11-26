package centaur.test.workflow

import java.nio.file.Path

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import centaur.test.metadata.WorkflowMetadata
import com.typesafe.config.{Config, ConfigFactory}
import configs.Result
import configs.syntax._
import cromwell.api.CromwellClient
import cromwell.api.model.WorkflowSingleSubmission
import common.validation.ErrorOr.ErrorOr

import scala.util.{Failure, Success, Try}

final case class Workflow private(testName: String,
                                  data: WorkflowData,
                                  metadata: Option[WorkflowMetadata],
                                  notInMetadata: List[String],
                                  directoryContentCounts: Option[DirectoryContentCountCheck],
                                  backends: BackendsRequirement) {
  def toWorkflowSubmission(refreshToken: Option[String]) = WorkflowSingleSubmission(
    wdl = data.wdl,
    workflowType = data.workflowType,
    workflowTypeVersion = data.workflowTypeVersion,
    inputsJson = data.inputs,
    options = CromwellClient.replaceJson(data.options, "refresh_token", refreshToken),
    labels = Option(data.labels),
    zippedImports = data.zippedImports)
}

object Workflow {

  def fromPath(path: Path): ErrorOr[Workflow] = {
    Try(ConfigFactory.parseFile(path.toFile)) match {
      case Success(c) => Workflow.fromConfig(c, path.getParent)
      case Failure(_) => invalidNel(s"Invalid test config: $path")
    }
  }

  def fromConfig(conf: Config, configPath: Path): ErrorOr[Workflow] = {
    conf.get[String]("name") match {
      case Result.Success(n) =>
        // If backend is provided, Centaur will only run this test if that backend is available on Cromwell
        val backendsRequirement = BackendsRequirement.fromConfig(conf.get[String]("backendsMode").map(_.toLowerCase).valueOrElse("all"), conf.get[List[String]]("backends").valueOrElse(List.empty[String]).map(_.toLowerCase))
        // If basePath is provided it'll be used as basis for finding other files, otherwise use the dir the config was in
        val basePath = conf.get[Path]("basePath") valueOrElse configPath
        val metadata: ErrorOr[Option[WorkflowMetadata]] = conf.get[Config]("metadata") match {
          case Result.Success(md) => WorkflowMetadata.fromConfig(md) map Option.apply
          case Result.Failure(_) => None.validNel
        }
        val absentMetadata = conf.get[List[String]]("absent-metadata-keys") match {
          case Result.Success(nmd) => nmd
          case Result.Failure(_) => List.empty
        }

        val directoryContentCheckValidation: ErrorOr[Option[DirectoryContentCountCheck]] = DirectoryContentCountCheck.forConfig(n, conf)
        val files = conf.get[Config]("files") match {
          case Result.Success(f) => WorkflowData.fromConfig(filesConfig = f, fullConfig = conf, basePath = basePath)
          case Result.Failure(_) => invalidNel(s"No 'files' block in $configPath")
        }

        (files, directoryContentCheckValidation, metadata) mapN {
          (f, d, m) => Workflow(n, f, m, absentMetadata, d, backendsRequirement)
        }

      case Result.Failure(_) => invalidNel(s"No name for: $configPath")
    }
  }
}
