package centaur.test.formulas

import cats.syntax.flatMap._
import cats.syntax.functor._
import centaur.test.Operations._
import centaur.test.Test.testMonad
import centaur.test.markers.CallMarker
import centaur.test.submit.{SubmitHttpResponse, SubmitResponse}
import centaur.test.workflow.Workflow
import centaur.test.{Operations, Test}
import centaur.{CentaurConfig, CromwellManager, CromwellTracker, ManagedCromwellServer}
import com.typesafe.config.ConfigFactory
import cromwell.api.model.{Aborted, Aborting, Failed, Running, SubmittedWorkflow, Succeeded, TerminalStatus, WorkflowMetadata}
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps

/**
  * A collection of test formulas which can be used, building upon operations by chaining them together via a
  * for comprehension. These assembled formulas can then be run by a client
  */
object TestFormulas {
  private def runWorkflowUntilTerminalStatus(workflow: Workflow, status: TerminalStatus): Test[SubmittedWorkflow] = {
    val workflowProgressTimeout = ConfigFactory.load().getOrElse("centaur.workflow-progress-timeout", 1 minute)
    for {
      s <- submitWorkflow(workflow)
      _ <- expectSomeProgress(s, workflow, Set(Running, status), workflowProgressTimeout)
      _ <- pollUntilStatus(s, workflow, status)
    } yield s
  }

  private def runSuccessfulWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = runWorkflowUntilTerminalStatus(workflow, Succeeded)
  private def runFailingWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = runWorkflowUntilTerminalStatus(workflow, Failed)

  def runSuccessfulWorkflowAndVerifyMetadata(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = for {
    submittedWorkflow <- runSuccessfulWorkflow(workflowDefinition)
    metadata <- validateMetadata(submittedWorkflow, workflowDefinition)
    _ = cromwellTracker.track(metadata)
    _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, metadata)
  } yield SubmitResponse(submittedWorkflow)

  def runFailingWorkflowAndVerifyMetadata(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = for {
    submittedWorkflow <- runFailingWorkflow(workflowDefinition)
    metadata <- validateMetadata(submittedWorkflow, workflowDefinition)
    _ = cromwellTracker.track(metadata)
    _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, metadata)
  } yield SubmitResponse(submittedWorkflow)

  def runWorkflowTwiceExpectingCaching(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    for {
      firstWF <- runSuccessfulWorkflow(workflowDefinition)
      secondWf <- runSuccessfulWorkflow(workflowDefinition.secondRun)
      _ <- printHashDifferential(firstWF, secondWf)
      metadata <- validateMetadata(secondWf, workflowDefinition, Option(firstWF.id.id))
      _ = cromwellTracker.track(metadata)
      _ <- validateNoCacheMisses(secondWf, metadata, workflowDefinition)
      _ <- validateDirectoryContentsCounts(workflowDefinition, secondWf, metadata)
    } yield SubmitResponse(secondWf)
  }

  def runWorkflowThriceExpectingCaching(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    for {
      firstWf <- runSuccessfulWorkflow(workflowDefinition)
      secondWf <- runSuccessfulWorkflow(workflowDefinition.secondRun)
      metadataTwo <- validateMetadata(secondWf, workflowDefinition, Option(firstWf.id.id))
      _ = cromwellTracker.track(metadataTwo)
      _ <- validateNoCacheHits(secondWf, metadataTwo, workflowDefinition)
      thirdWf <- runSuccessfulWorkflow(workflowDefinition.thirdRun)
      _ <- printHashDifferential(secondWf, thirdWf)
      metadataThree <- validateMetadata(thirdWf, workflowDefinition, Option(secondWf.id.id))
      _ <- validateNoCacheMisses(thirdWf, metadataThree, workflowDefinition)
      _ <- validateDirectoryContentsCounts(workflowDefinition, thirdWf, metadataThree)
    } yield SubmitResponse(thirdWf)
  }

  def runWorkflowTwiceExpectingNoCaching(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    for {
      _ <- runSuccessfulWorkflow(workflowDefinition) // Build caches
      testWf <- runSuccessfulWorkflow(workflowDefinition.secondRun)
      metadata <- validateMetadata(testWf, workflowDefinition)
      _ = cromwellTracker.track(metadata)
      _ <- validateNoCacheHits(testWf, metadata, workflowDefinition)
      _ <- validateDirectoryContentsCounts(workflowDefinition, testWf, metadata)
    } yield SubmitResponse(testWf)
  }

  def runFailingWorkflowTwiceExpectingNoCaching(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    for {
      _ <- runFailingWorkflow(workflowDefinition) // Build caches
      testWf <- runFailingWorkflow(workflowDefinition)
      metadata <- validateMetadata(testWf, workflowDefinition)
      _ = cromwellTracker.track(metadata)
      _ <- validateNoCacheHits(testWf, metadata, workflowDefinition)
      _ <- validateDirectoryContentsCounts(workflowDefinition, testWf, metadata)
    } yield SubmitResponse(testWf)
  }

  private def cromwellRestart(workflowDefinition: Workflow,
                              callMarker: CallMarker,
                              testRecover: Boolean,
                              finalStatus: TerminalStatus)(
                              implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    CentaurConfig.runMode match {
      case ManagedCromwellServer(_, postRestart, withRestart) if withRestart =>
        for {
          submittedWorkflow <- submitWorkflow(workflowDefinition)
          jobId <- pollUntilCallIsRunning(workflowDefinition, submittedWorkflow, callMarker.callKey)
          _ = CromwellManager.stopCromwell(s"Scheduled restart from ${workflowDefinition.testName}")
          _ = CromwellManager.startCromwell(postRestart)
          _ <- expectSomeProgress(submittedWorkflow, workflowDefinition, Set(Running, finalStatus), 1.minute)
          _ <- pollUntilStatus(submittedWorkflow, workflowDefinition, finalStatus)
          metadata <- validateMetadata(submittedWorkflow, workflowDefinition)
          _ = cromwellTracker.track(metadata)
          _ <- if (testRecover) {
            validateRecovered(workflowDefinition, submittedWorkflow, metadata, callMarker.callKey, jobId)
          }
          else {
            Test.successful(())
          }
          _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, metadata)
        } yield SubmitResponse(submittedWorkflow)
      case _ if finalStatus == Succeeded => runSuccessfulWorkflowAndVerifyMetadata(workflowDefinition)
      case _ if finalStatus == Failed => runFailingWorkflowAndVerifyMetadata(workflowDefinition)
      case _ => Test.invalidTestDefinition("This test can only run successful or failed workflow", workflowDefinition)
    }
  }

  def instantAbort(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = for {
    submittedWorkflow <- submitWorkflow(workflowDefinition)
    _ <- abortWorkflow(submittedWorkflow)
    _ <- expectSomeProgress(submittedWorkflow, workflowDefinition, Set(Running, Aborting, Aborted), 1.minute)
    _ <- pollUntilStatus(submittedWorkflow, workflowDefinition, Aborted)
    metadata <- validateMetadata(submittedWorkflow, workflowDefinition)
    _ = cromwellTracker.track(metadata)
    _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, metadata)
  } yield SubmitResponse(submittedWorkflow)

  def scheduledAbort(workflowDefinition: Workflow, callMarker: CallMarker, restart: Boolean)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    def withRestart() = CentaurConfig.runMode match {
      case ManagedCromwellServer(_, postRestart, withRestart) if withRestart =>
        CromwellManager.stopCromwell(s"Scheduled restart from ${workflowDefinition.testName}")
        CromwellManager.startCromwell(postRestart)
      case _ =>
    }
    
    for {
      submittedWorkflow <- submitWorkflow(workflowDefinition)
      jobId <- pollUntilCallIsRunning(workflowDefinition, submittedWorkflow, callMarker.callKey)
      // The Cromwell call status could be running but the backend job might not have started yet, give it some time
      _ <- waitFor(30.seconds)
      _ <- abortWorkflow(submittedWorkflow)
      _ = if(restart) withRestart()
      _ <- expectSomeProgress(submittedWorkflow, workflowDefinition, Set(Running, Aborting, Aborted), 1.minute)
      _ <- pollUntilStatus(submittedWorkflow, workflowDefinition, Aborted)
      _ <- validatePAPIAborted(workflowDefinition, submittedWorkflow, jobId)
      // Wait a little to make sure that if the abort didn't work and calls start running we see them in the metadata
      _ <- waitFor(30.seconds)
      metadata <- validateMetadata(submittedWorkflow, workflowDefinition)
      _ = cromwellTracker.track(metadata)
      _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, metadata)
    } yield SubmitResponse(submittedWorkflow)
  }

  def workflowRestart(workflowDefinition: Workflow,
                             callMarker: CallMarker,
                             recover: Boolean,
                             finalStatus: TerminalStatus)(
                             implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    cromwellRestart(workflowDefinition, callMarker, testRecover = recover, finalStatus = finalStatus)
  }

  def submitInvalidWorkflow(workflow: Workflow, expectedSubmitResponse: SubmitHttpResponse): Test[SubmitResponse] = {
    for {
      actualSubmitResponse <- Operations.submitInvalidWorkflow(workflow)
      _ <- validateSubmitFailure(workflow, expectedSubmitResponse, actualSubmitResponse)
    } yield actualSubmitResponse
  }

  def papiUpgrade(workflowDefinition: Workflow, callMarker: CallMarker)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    CentaurConfig.runMode match {
      case ManagedCromwellServer(_, postRestart, withRestart) if withRestart =>
        for {
          first <- submitWorkflow(workflowDefinition)
          _ <- pollUntilCallIsRunning(workflowDefinition, first, callMarker.callKey)
          _ = CromwellManager.stopCromwell(s"Scheduled restart from ${workflowDefinition.testName}")
          _ = CromwellManager.startCromwell(postRestart)
          _ <- expectSomeProgress(first, workflowDefinition, Set(Running, Succeeded), 1.minute)
          _ <- pollUntilStatus(first, workflowDefinition, Succeeded)
          second <- runSuccessfulWorkflow(workflowDefinition.secondRun) // Same WDL and config but a "backend" runtime option targeting PAPI v2.
          _ <- printHashDifferential(first, second)
          metadata <- validateMetadata(second, workflowDefinition, Option(first.id.id))
          _ = cromwellTracker.track(metadata)
          _ <- validateNoCacheMisses(second, metadata, workflowDefinition)
          _ <- validateDirectoryContentsCounts(workflowDefinition, second, metadata)
        } yield SubmitResponse(second)
      case _ => Test.invalidTestDefinition("Configuration not supported by PapiUpgradeTest", workflowDefinition)
    }
  }

  implicit class EnhancedCromwellTracker(val tracker: Option[CromwellTracker]) extends AnyVal {
    def track(metadata: WorkflowMetadata): Unit = tracker foreach { _.track(metadata) }
  }
}
