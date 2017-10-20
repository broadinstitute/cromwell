package centaur.test.standard

import java.nio.file.Path

import cats.data.Validated._
import cats.implicits._
import centaur.test._
import centaur.test.formulas.TestFormulas
import centaur.test.standard.CentaurTestFormat._
import centaur.test.submit.SubmitResponse
import centaur.test.workflow.{AllBackendsRequired, AnyBackendRequired, OnlyBackendsAllowed, Workflow}
import com.typesafe.config.{Config, ConfigFactory}
import lenthall.validation.ErrorOr.ErrorOr

import scala.util.{Failure, Success, Try}

case class CentaurTestCase(workflow: Workflow,
                           testFormat: CentaurTestFormat,
                           testOptions: TestOptions,
                           submitResponseOption: Option[SubmitResponse]) {
  def testFunction: Test[Unit] = this.testFormat match {
    case WorkflowSuccessTest => TestFormulas.runSuccessfulWorkflowAndVerifyMetadata(workflow)
    case WorkflowFailureTest => TestFormulas.runFailingWorkflowAndVerifyMetadata(workflow)
    case RunTwiceExpectingCallCachingTest => TestFormulas.runWorkflowTwiceExpectingCaching(workflow)
    case RunTwiceExpectingNoCallCachingTest => TestFormulas.runWorkflowTwiceExpectingNoCaching(workflow)
    case RunFailingTwiceExpectingNoCallCachingTest => TestFormulas.runFailingWorkflowTwiceExpectingNoCaching(workflow)
    case SubmitFailureTest => TestFormulas.submitInvalidWorkflow(workflow, submitResponseOption.get)
    case InstantAbort => TestFormulas.instantAbort(workflow)
    case CromwellRestartWithRecover(callMarker)=> TestFormulas.cromwellRestartWithRecover(workflow, callMarker)
    case CromwellRestartWithoutRecover(callMarker) => TestFormulas.cromwellRestartWithoutRecover(workflow, callMarker)
    case ScheduledAbort(callMarker) => TestFormulas.scheduledAbort(workflow, callMarker, restart = false)
    case ScheduledAbortWithRestart(callMarker) => TestFormulas.scheduledAbort(workflow, callMarker, restart = true)
    case other => Test.failed(new Exception(s"Invalid test format $other"))
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
      case Success(c) =>
        CentaurTestCase.fromConfig(c, path.getParent) match {
          case Valid(testCase) => validateTestCase(testCase)
          case invalid: Invalid[_] => invalid
        }
      case Failure(f) => invalidNel(s"Invalid test config: $path (${f.getMessage})")
    }
  }

  def fromConfig(conf: Config, configPath: Path): ErrorOr[CentaurTestCase] = {
    val workflow = Workflow.fromConfig(conf, configPath)
    val format = CentaurTestFormat.fromConfig(conf).toValidated
    val options = TestOptions.fromConfig(conf)
    val submit = SubmitResponse.fromConfig(conf)
    workflow |@| format |@| options |@| submit map {
      CentaurTestCase(_, _, _, _)
    }
  }

  private def validateTestCase(testCase: CentaurTestCase): ErrorOr[CentaurTestCase] = {
    testCase.testFormat match {
      case SubmitFailureTest => validateSubmitFailure(testCase.workflow, testCase.submitResponseOption).map(_ => testCase)
      case _ => Valid(testCase)
    }
  }

  private def validateSubmitFailure(workflow: Workflow, submitResponseOption: Option[SubmitResponse]): ErrorOr[Unit] = {
    submitResponseOption match {
      case None => invalidNel("No submit stanza included in test config")
      case Some(_) => Valid(())
    }
  }
}
