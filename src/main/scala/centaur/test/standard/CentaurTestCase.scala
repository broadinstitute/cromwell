package centaur.test.standard

import java.nio.file.Path

import cats.data.Validated._
import cats.Apply
import centaur.test._
import centaur.test.formulas.TestFormulas
import centaur.test.standard.CentaurTestFormat._
import centaur.test.workflow.{AllBackendsRequired, AnyBackendRequired, OnlyBackendsAllowed, Workflow}
import com.typesafe.config.{Config, ConfigFactory}

import scala.util.{Failure, Success, Try}

case class CentaurTestCase(workflow: Workflow, testFormat: CentaurTestFormat, testOptions: TestOptions) {
  def testFunction: Test[Unit] = this.testFormat match {
    case WorkflowSuccessTest => TestFormulas.runSuccessfulWorkflowAndVerifyMetadata(workflow)
    case WorkflowFailureTest => TestFormulas.runFailingWorkflowAndVerifyMetadata(workflow)
    case RunTwiceExpectingCallCachingTest => TestFormulas.runWorkflowTwiceExpectingCaching(workflow)
    case RunTwiceExpectingNoCallCachingTest => TestFormulas.runWorkflowTwiceExpectingNoCaching(workflow)
    case RunFailingTwiceExpectingNoCallCachingTest => TestFormulas.runFailingWorkflowTwiceExpectingNoCaching(workflow)
    case CromwellRestartWithRecover => TestFormulas.cromwellRestartWithRecover(workflow)
    case CromwellRestartWithoutRecover => TestFormulas.cromwellRestartWithoutRecover(workflow)
  }

  def isIgnored(supportedBackends: List[String]): Boolean = {
    val backendSupported = workflow.backends match {
      case AllBackendsRequired(allBackends) => allBackends forall supportedBackends.contains
      case AnyBackendRequired(anyBackend) => anyBackend exists supportedBackends.contains
      case OnlyBackendsAllowed(onlyBackends) => onlyBackends.toSet == supportedBackends.toSet
    }

    testOptions.ignore || !backendSupported
  }
}

object CentaurTestCase {
  def fromPath(path: Path): ErrorOr[CentaurTestCase] = {
    Try(ConfigFactory.parseFile(path.toFile)) match {
      case Success(c) => CentaurTestCase.fromConfig(c, path.getParent)
      case Failure(f) => invalidNel(s"Invalid test config: $path (${f.getMessage})")
    }
  }

  def fromConfig(conf: Config, configPath: Path): ErrorOr[CentaurTestCase] = {
    val workflow = Workflow.fromConfig(conf, configPath)
    val format = CentaurTestFormat.fromConfig(conf)
    val options = TestOptions.fromConfig(conf)
    Apply[ErrorOr].map3(workflow, format, options)((w, f, o) => CentaurTestCase(w, f, o))
  }
}
