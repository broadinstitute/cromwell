package centaur.test.standard

import java.nio.file.Path

import cats.data.Validated._
import cats.Apply
import cats.std.list._
import centaur.test._
import centaur.test.formulas.TestFormulas
import centaur.test.standard.StandardTestFormat.{SubmissionFailureTest, WorkflowFailureTest, WorkflowSuccessTest}
import centaur.test.workflow.Workflow
import centaur.test.workflow.Workflow.{WorkflowWithMetadata, WorkflowWithoutMetadata}
import com.typesafe.config.{Config, ConfigFactory}

import scala.util.{Failure, Success, Try}

case class StandardTestCase(workflow: Workflow, testFormat: StandardTestFormat) {
  def testFunction = this.testFormat match {
    case WorkflowSuccessTest => successfulTestFunction
    case WorkflowFailureTest => TestFormulas.runFailingWorkflow _
    case SubmissionFailureTest => TestFormulas.runSubmissionFailureWorkflow _
  }

  private def successfulTestFunction = this.workflow match {
    case _: WorkflowWithoutMetadata => TestFormulas.runSuccessfulWorkflow _
    case _: WorkflowWithMetadata => TestFormulas.runSuccessfulWorkflowAndVerifyMetadata _
  }
}

object StandardTestCase {
  def fromPath(path: Path): ErrorOr[StandardTestCase] = {
    Try(ConfigFactory.parseFile(path.toFile)) match {
      case Success(c) => StandardTestCase.fromConfig(c, path.getParent)
      case Failure(f) => invalidNel(s"Invalid test config: $path")
    }
  }

  def fromConfig(conf: Config, configPath: Path): ErrorOr[StandardTestCase] = {
    val workflow = Workflow.fromConfig(conf, configPath)
    val format = StandardTestFormat.fromConfig(conf)
    Apply[ErrorOr].map2(workflow, format)((w, f) => StandardTestCase(w, f))
  }
}
