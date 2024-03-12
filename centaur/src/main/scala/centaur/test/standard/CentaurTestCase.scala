package centaur.test.standard

import better.files._
import cats.data.Validated._
import cats.syntax.all._
import centaur.CromwellTracker
import centaur.test._
import centaur.test.formulas.TestFormulas
import centaur.test.standard.CentaurTestFormat._
import centaur.test.submit.{SubmitHttpResponse, SubmitResponse}
import centaur.test.workflow._
import com.typesafe.config.{Config, ConfigFactory}
import common.validation.ErrorOr._
import cromwell.api.model.{Failed, Succeeded}

import scala.util.{Failure, Success, Try}

case class CentaurTestCase(workflow: Workflow,
                           testFormat: CentaurTestFormat,
                           testOptions: TestOptions,
                           submittedWorkflowTracker: SubmittedWorkflowTracker,
                           submitResponseOption: Option[SubmitHttpResponse]
)(implicit cromwellTracker: Option[CromwellTracker]) {

  def testFunction: Test[SubmitResponse] = this.testFormat match {
    case WorkflowSuccessTest => TestFormulas.runSuccessfulWorkflowAndVerifyMetadata(workflow)
    case WorkflowSuccessAndTimedOutputsTest => TestFormulas.runSuccessfulWorkflowAndVerifyTimeAndOutputs(workflow)
    case WorkflowFailureTest => TestFormulas.runFailingWorkflowAndVerifyMetadata(workflow)
    case RunTwiceExpectingCallCachingTest => TestFormulas.runWorkflowTwiceExpectingCaching(workflow)
    case RunThriceExpectingCallCachingTest => TestFormulas.runWorkflowThriceExpectingCaching(workflow)
    case RunTwiceExpectingNoCallCachingTest => TestFormulas.runWorkflowTwiceExpectingNoCaching(workflow)
    case RunFailingTwiceExpectingNoCallCachingTest => TestFormulas.runFailingWorkflowTwiceExpectingNoCaching(workflow)
    case SubmitFailureTest => TestFormulas.submitInvalidWorkflow(workflow, submitResponseOption.get)
    case InstantAbort => TestFormulas.instantAbort(workflow)
    case CromwellRestartWithRecover(callMarker) =>
      TestFormulas.workflowRestart(workflow, callMarker, recover = true, finalStatus = Succeeded)
    case WorkflowFailureRestartWithRecover(callMarker) =>
      TestFormulas.workflowRestart(workflow, callMarker, recover = true, finalStatus = Failed)
    case WorkflowFailureRestartWithoutRecover(callMarker) =>
      TestFormulas.workflowRestart(workflow, callMarker, recover = false, finalStatus = Failed)
    case CromwellRestartWithoutRecover(callMarker) =>
      TestFormulas.workflowRestart(workflow, callMarker, recover = false, finalStatus = Succeeded)
    case ScheduledAbort(callMarker) => TestFormulas.scheduledAbort(workflow, callMarker, restart = false)
    case ScheduledAbortWithRestart(callMarker) => TestFormulas.scheduledAbort(workflow, callMarker, restart = true)
    case PapiUpgradeTest(callMarker) => TestFormulas.papiUpgrade(workflow, callMarker)
    case other => Test.invalidTestDefinition(s"Invalid test format $other", workflow)
  }

  def isIgnored(supportedBackends: List[String]): Boolean = {
    val backendSupported = workflow.backends match {
      case AllBackendsRequired(testBackends) =>
        // Test will run on servers that support all of the test's backends (or more) (default)
        testBackends forall supportedBackends.contains
      case AnyBackendRequired(testBackends) =>
        // Test will run on servers that support at least one of the test's backends (or more)
        testBackends exists supportedBackends.contains
      case OnlyBackendsAllowed(testBackends) =>
        // Test will run on servers that only support backends the test specifies (or fewer)
        supportedBackends forall testBackends.contains
    }

    testOptions.ignore || !backendSupported
  }

  def containsTag(tag: String): Boolean = testOptions.tags.contains(tag)

  def name: String = s"${testFormat.testSpecString} ${workflow.testName}"
}

object CentaurTestCase {
  def fromFile(cromwellTracker: Option[CromwellTracker])(file: File): ErrorOr[CentaurTestCase] =
    Try(ConfigFactory.parseFile(file.toJava).resolve()) match {
      case Success(c) =>
        CentaurTestCase.fromConfig(c, file.parent, cromwellTracker) flatMap validateTestCase leftMap {
          s"Error in test file '$file'." :: _
        }
      case Failure(f) => invalidNel(s"Invalid test config: $file (${f.getMessage})")
    }

  def fromConfig(conf: Config, configFile: File, cromwellTracker: Option[CromwellTracker]): ErrorOr[CentaurTestCase] = {
    val submittedWorkflowTracker = new SubmittedWorkflowTracker()
    val workflow = Workflow.fromConfig(conf, configFile, submittedWorkflowTracker)
    val format: ErrorOr[CentaurTestFormat] = CentaurTestFormat.fromConfig(conf).toValidated
    val options = TestOptions.fromConfig(conf)
    val submit = SubmitHttpResponse.fromConfig(conf)
    (workflow, format, options, submit) mapN {
      CentaurTestCase(_, _, _, submittedWorkflowTracker, _)(cromwellTracker)
    }
  }

  private def validateTestCase(testCase: CentaurTestCase): ErrorOr[CentaurTestCase] =
    testCase.testFormat match {
      case SubmitFailureTest =>
        validateSubmitFailure(testCase.workflow, testCase.submitResponseOption).map(_ => testCase)
      case _ => Valid(testCase)
    }

  private def validateSubmitFailure(workflow: Workflow,
                                    submitResponseOption: Option[SubmitHttpResponse]
  ): ErrorOr[SubmitResponse] =
    submitResponseOption match {
      case None => invalidNel("No submit stanza included in test config")
      case Some(response) => Valid(response)
    }
}
