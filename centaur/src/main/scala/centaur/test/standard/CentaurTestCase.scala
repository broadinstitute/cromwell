package centaur.test.standard

import better.files._
import cats.data.Validated._
import cats.effect.IO
import cats.syntax.all._
import centaur.CromwellTracker
import centaur.test._
import centaur.test.formulas.TestFormulas
import centaur.test.standard.CentaurTestFormat._
import centaur.test.submit.{SubmitHttpResponse, SubmitResponse}
import centaur.test.workflow.{AllBackendsRequired, AnyBackendRequired, OnlyBackendsAllowed, Workflow}
import com.typesafe.config.{Config, ConfigFactory}
import common.validation.ErrorOr._
import cromwell.api.model.{Failed, SubmittedWorkflow, Succeeded, WorkflowId}

import scala.util.{Failure, Success, Try}

case class CentaurTestCase(workflow: Workflow,
                           testFormat: CentaurTestFormat,
                           testOptions: TestOptions,
                           submitResponseOption: Option[SubmitHttpResponse])(
                           implicit cromwellTracker: Option[CromwellTracker]) {

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
    case CromwellRestartWithRecover(callMarker)=> TestFormulas.workflowRestart(workflow, callMarker, recover = true, finalStatus = Succeeded)
    case WorkflowFailureRestartWithRecover(callMarker)=> TestFormulas.workflowRestart(workflow, callMarker, recover = true, finalStatus = Failed)
    case WorkflowFailureRestartWithoutRecover(callMarker)=> TestFormulas.workflowRestart(workflow, callMarker, recover = false, finalStatus = Failed)
    case CromwellRestartWithoutRecover(callMarker) => TestFormulas.workflowRestart(workflow, callMarker, recover = false, finalStatus = Succeeded)
    case ScheduledAbort(callMarker) => TestFormulas.scheduledAbort(workflow, callMarker, restart = false)
    case ScheduledAbortWithRestart(callMarker) => TestFormulas.scheduledAbort(workflow, callMarker, restart = true)
    case PapiUpgradeTest(callMarker) => TestFormulas.papiUpgrade(workflow, callMarker)
    case other => Test.invalidTestDefinition(s"Invalid test format $other", workflow)
  }

  def isIgnored(supportedBackends: List[String]): Boolean = {
    val backendSupported = workflow.backends match {
      case AllBackendsRequired(allBackends) => allBackends forall supportedBackends.contains
      case AnyBackendRequired(anyBackend) => anyBackend exists supportedBackends.contains
      case OnlyBackendsAllowed(onlyBackends) => supportedBackends forall onlyBackends.contains
    }

    testOptions.ignore || !backendSupported
  }

  def containsTag(tag: String): Boolean = testOptions.tags.contains(tag)

  def name: String = s"${testFormat.testSpecString} ${workflow.testName}"

  private var submittedWorkflowIds: List[WorkflowId] = List.empty

  /**
   * Run the specified cleanup function on the submitted workflow IDs tracked by this `CentaurTestCase`, clearing out
   * the list of submitted workflow IDs afterward.
   */
  def cleanUpBeforeRetry(cleanUpFunction: WorkflowId => IO[Unit]): IO[Unit] = for {
    _ <- submittedWorkflowIds.traverse(cleanUpFunction)
    _ = submittedWorkflowIds = List.empty
  } yield ()

  /**
   * Add a `SubmittedWorkflow` to the list of `SubmittedWorkflow`s to clean up should this `CentaurTestCase` require a
   * retry. Prevents unwanted cache hits from partially successful attempts when retrying a call caching test case.
   */
  def addSubmittedWorkflow(submittedWorkflow: SubmittedWorkflow): Unit = {
    submittedWorkflowIds = submittedWorkflow.id :: submittedWorkflowIds
  }
}

object CentaurTestCase {
  def fromFile(cromwellTracker: Option[CromwellTracker])(file: File): ErrorOr[CentaurTestCase] = {
    Try(ConfigFactory.parseFile(file.toJava).resolve()) match {
      case Success(c) =>
        CentaurTestCase.fromConfig(c, file.parent, cromwellTracker) flatMap validateTestCase leftMap { s"Error in test file '$file'." :: _ }
      case Failure(f) => invalidNel(s"Invalid test config: $file (${f.getMessage})")
    }
  }

  def fromConfig(conf: Config, configFile: File, cromwellTracker: Option[CromwellTracker]): ErrorOr[CentaurTestCase] = {
    val workflow = Workflow.fromConfig(conf, configFile)
    val format: ErrorOr[CentaurTestFormat] = CentaurTestFormat.fromConfig(conf).toValidated
    val options = TestOptions.fromConfig(conf)
    val submit = SubmitHttpResponse.fromConfig(conf)
    for {
      testCase <- (workflow, format, options, submit) mapN { CentaurTestCase(_, _, _, _)(cromwellTracker) }
      _ = testCase.workflow.setTestCase(testCase)
    } yield testCase
  }

  private def validateTestCase(testCase: CentaurTestCase): ErrorOr[CentaurTestCase] = {
    testCase.testFormat match {
      case SubmitFailureTest => validateSubmitFailure(testCase.workflow, testCase.submitResponseOption).map(_ => testCase)
      case _ => Valid(testCase)
    }
  }

  private def validateSubmitFailure(workflow: Workflow,
                                    submitResponseOption: Option[SubmitHttpResponse]): ErrorOr[SubmitResponse] = {
    submitResponseOption match {
      case None => invalidNel("No submit stanza included in test config")
      case Some(response) => Valid(response)
    }
  }
}
