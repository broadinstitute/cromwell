package centaur.test.workflow

import java.nio.file.Path

import cats.Apply
import cats.data.Validated._
import centaur.test._
import centaur.test.metadata.WorkflowMetadata
import com.typesafe.config.{Config, ConfigFactory}
import configs.Result
import configs.syntax._

import scala.util.{Failure, Success, Try}

sealed abstract class Workflow {
  def testName: String
  def data: WorkflowData
  def backends: List[String]
}

object Workflow {
  final case class WorkflowWithMetadata(testName: String, data: WorkflowData, metadata: WorkflowMetadata, backends: List[String]) extends Workflow
  final case class WorkflowWithoutMetadata(testName: String, data: WorkflowData, backends: List[String]) extends Workflow

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
        val backends = conf.get[List[String]]("backends").valueOrElse(List.empty[String]).map(_.toLowerCase)
        // If basePath is provided it'll be used as basis for finding other files, otherwise use the dir the config was in
        val basePath = conf.get[Path]("basePath") valueOrElse configPath
        val metadata = conf.get[Config]("metadata")
        val files = conf.get[Config]("files") match {
          case Result.Success(c) => WorkflowData.fromConfig(c, basePath)
          case Result.Failure(_) => invalidNel(s"No 'files' block in $configPath")
        }

        metadata match {
          case Result.Success(c) => Apply[ErrorOr].map2(files, WorkflowMetadata.fromConfig(c))((f, m) => WorkflowWithMetadata(n, f, m, backends))
          case Result.Failure(_) => Apply[ErrorOr].map(files)(f => WorkflowWithoutMetadata(n, f, backends))
        }

      case Result.Failure(_) => invalidNel(s"No name for: $configPath")
    }
  }
}
