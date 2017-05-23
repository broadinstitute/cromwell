package centaur.test.standard

import java.nio.file.Path

import cats.data.Validated._
import cats.Apply
import centaur.test._
import centaur.test.formulas.TestFormulas
import centaur.test.standard.StandardTestFormat.{WorkflowFailureTest, WorkflowSuccessTest}
import centaur.test.workflow.Workflow
import centaur.test.workflow.Workflow.{WorkflowWithMetadata, WorkflowWithoutMetadata}
import com.typesafe.config.{Config, ConfigFactory}

import scala.util.{Failure, Success, Try}

case class StandardTestCase(workflow: Workflow, testFormat: StandardTestFormat, testOptions: TestOptions) {
  def testFunction = this.testFormat match {
    case WorkflowSuccessTest => successfulTestFunction
    case WorkflowFailureTest => failureTestFunction
  }

  def isIgnored(supportedBackends: List[String]): Boolean = {
    val backendSupported = workflow.backends forall supportedBackends.contains
    testOptions.ignore || !backendSupported
  }

  private def successfulTestFunction = this.workflow match {
    case _: WorkflowWithoutMetadata => TestFormulas.runSuccessfulWorkflow _
    case _: WorkflowWithMetadata => TestFormulas.runSuccessfulWorkflowAndVerifyMetadata _
  }

  private def failureTestFunction = this.workflow match {
    case _: WorkflowWithoutMetadata => TestFormulas.runFailingWorkflow _
    case _: WorkflowWithMetadata => TestFormulas.runFailingWorkflowAndVerifyMetadata _
  }
}

object StandardTestCase {
  def fromPath(path: Path): ErrorOr[StandardTestCase] = {
    Try(ConfigFactory.parseFile(path.toFile)) match {
      case Success(c) => StandardTestCase.fromConfig(c, path.getParent)
      case Failure(f) => invalidNel(s"Invalid test config: $path due to failure: $f")
    }
  }

  def fromConfig(conf: Config, configPath: Path): ErrorOr[StandardTestCase] = {
    val workflow = Workflow.fromConfig(conf, configPath)
    val format = StandardTestFormat.fromConfig(conf)
    val options = TestOptions.fromConfig(conf)
    Apply[ErrorOr].map3(workflow, format, options)((w, f, o) => StandardTestCase(w, f, o))
  }
}
