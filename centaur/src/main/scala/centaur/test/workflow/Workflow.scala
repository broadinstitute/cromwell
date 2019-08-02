package centaur.test.workflow

import java.nio.file.Path

import common.validation.Validation._
import better.files._
import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import centaur.test.metadata.WorkflowFlatMetadata
import com.typesafe.config.{Config, ConfigFactory}
import common.validation.ErrorOr.ErrorOr
import configs.Result
import configs.syntax._
import cromwell.api.CromwellClient
import cromwell.api.model.WorkflowSingleSubmission

import scala.util.{Failure, Success, Try}

final case class Workflow private(testName: String,
                                  data: WorkflowData,
                                  metadata: Option[WorkflowFlatMetadata],
                                  notInMetadata: List[String],
                                  directoryContentCounts: Option[DirectoryContentCountCheck],
                                  backends: BackendsRequirement,
                                  retryTestFailures: Boolean,
                                  allowOtherOutputs: Boolean) {
  def toWorkflowSubmission(refreshToken: Option[String]) = WorkflowSingleSubmission(
    workflowSource = data.workflowContent,
    workflowUrl = data.workflowUrl,
    workflowRoot = data.workflowRoot,
    workflowType = data.workflowType,
    workflowTypeVersion = data.workflowTypeVersion,
    inputsJson = data.inputs.map(_.unsafeRunSync()),
    options = CromwellClient.replaceJson(data.options.map(_.unsafeRunSync()), "refresh_token", refreshToken),
    labels = Option(data.labels),
    zippedImports = data.zippedImports)

  def secondRun: Workflow = {
    copy(data = data.copy(options = data.secondOptions))
  }

  def thirdRun: Workflow = {
    copy(data = data.copy(options = data.thirdOptions))
  }
}

object Workflow {

  def fromPath(path: Path): ErrorOr[Workflow] = {
    Try(ConfigFactory.parseFile(path.toFile).resolve()) match {
      case Success(c) => Workflow.fromConfig(c, path.getParent)
      case Failure(_) => invalidNel(s"Invalid test config: $path")
    }
  }

  def fromConfig(conf: Config, configFile: File): ErrorOr[Workflow] = {
    conf.get[String]("name") match {
      case Result.Success(n) =>
        // If backend is provided, Centaur will only run this test if that backend is available on Cromwell
        val backendsRequirement = BackendsRequirement.fromConfig(conf.get[String]("backendsMode").map(_.toLowerCase).valueOrElse("all"), conf.get[List[String]]("backends").valueOrElse(List.empty[String]).map(_.toLowerCase))
        // If basePath is provided it'll be used as basis for finding other files, otherwise use the dir the config was in
        val basePath = conf.get[Option[Path]]("basePath") valueOrElse None map (File(_)) getOrElse configFile
        val metadata: ErrorOr[Option[WorkflowFlatMetadata]] = conf.get[Config]("metadata") match {
          case Result.Success(md) => WorkflowFlatMetadata.fromConfig(md) map Option.apply
          case Result.Failure(_) => None.validNel
        }
        val absentMetadata = conf.get[List[String]]("absent-metadata-keys") match {
          case Result.Success(nmd) => nmd
          case Result.Failure(_) => List.empty
        }

        val directoryContentCheckValidation: ErrorOr[Option[DirectoryContentCountCheck]] = DirectoryContentCountCheck.forConfig(n, conf)
        val files = conf.get[Config]("files") match {
          case Result.Success(f) => WorkflowData.fromConfig(filesConfig = f, fullConfig = conf, basePath = basePath)
          case Result.Failure(_) => invalidNel(s"No 'files' block in $configFile")
        }
        val retryTestFailuresErrorOr = validate(conf.get[Boolean]("retryTestFailures").valueOrElse(true))
        val allowOtherOutputs: Boolean = conf.get[Boolean]("allowOtherOutputs") match {
          case Result.Success(allow) => allow
          case Result.Failure(_) => true
        }

        (files, directoryContentCheckValidation, metadata, retryTestFailuresErrorOr) mapN {
          (f, d, m, retryTestFailures) => Workflow(n, f, m, absentMetadata, d, backendsRequirement, retryTestFailures, allowOtherOutputs)
        }

      case Result.Failure(_) => invalidNel(s"No name for: $configFile")
    }
  }
}
