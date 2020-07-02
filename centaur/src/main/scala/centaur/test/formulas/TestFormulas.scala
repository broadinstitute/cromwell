package centaur.test.formulas

import cats.syntax.flatMap._
import cats.syntax.functor._
import centaur.api.CentaurCromwellClient
import centaur.test.Operations._
import centaur.test.Test._
import centaur.test.markers.CallMarker
import centaur.test.submit.{SubmitHttpResponse, SubmitResponse}
import centaur.test.workflow.Workflow
import centaur.test.{Operations, Test}
import centaur.{CentaurConfig, CromwellManager, CromwellTracker, ManagedCromwellServer}
import com.typesafe.config.ConfigFactory
import com.typesafe.scalalogging.StrictLogging
import cromwell.api.model.{Aborted, Aborting, Failed, Running, SubmittedWorkflow, Succeeded, TerminalStatus, WorkflowMetadata}
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps
import centaur.test.metadata.WorkflowFlatMetadata._
import spray.json.JsString
/**
  * A collection of test formulas which can be used, building upon operations by chaining them together via a
  * for comprehension. These assembled formulas can then be run by a client
  */
object TestFormulas extends StrictLogging {
  val workflowProgressTimeout = ConfigFactory.load().getOrElse("centaur.workflow-progress-timeout", 1 minute)
  logger.info(s"Running with a workflow progress timeout of $workflowProgressTimeout")

  private def runWorkflowUntilTerminalStatus(workflow: Workflow, status: TerminalStatus): Test[SubmittedWorkflow] = {
    for {
      _ <- checkVersion()
      s <- submitWorkflow(workflow)
      _ <- expectSomeProgress(s, workflow, Set(Running, status), workflowProgressTimeout)
      _ <- pollUntilStatus(s, workflow, status)
    } yield s
  }

  private def runSuccessfulWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = runWorkflowUntilTerminalStatus(workflow, Succeeded)
  private def runFailingWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = runWorkflowUntilTerminalStatus(workflow, Failed)

  def runSuccessfulWorkflowAndVerifyMetadata(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = for {
    _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
    submittedWorkflow <- runSuccessfulWorkflow(workflowDefinition)
    labelsLikelyBeforeArchival = CentaurCromwellClient.labels(submittedWorkflow)
    unarchivedNonSubworkflowMetadata <- fetchAndValidateNonSubworkflowMetadata(submittedWorkflow, workflowDefinition, validateArchived = Option(false))
    unarchivedFullMetadata <- fetchMetadata(submittedWorkflow, expandSubworkflows = true, requestArchivedMetadata = Option(false)).asTest
    unarchivedJobManagerStyleMetadata <- fetchAndValidateJobManagerStyleMetadata(submittedWorkflow, workflowDefinition, prefetchedOriginalNonSubWorkflowMetadata = Option(unarchivedNonSubworkflowMetadata.value), validateArchived = Option(false))
    notArchivedFlatMetadata = unarchivedFullMetadata.asFlat
    workflowRoot = notArchivedFlatMetadata.value.get("workflowRoot").collectFirst { case JsString(r) => r } getOrElse "No Workflow Root"
    unarchivedOutputs <- fetchAndValidateOutputs(submittedWorkflow, workflowDefinition, workflowRoot, validateArchived = Option(false))
    _ <- fetchAndValidateLabels(submittedWorkflow, workflowDefinition, workflowRoot)
    _ <- validateLogs(unarchivedNonSubworkflowMetadata, submittedWorkflow, workflowDefinition, validateArchived = Option(false))
    _ = cromwellTracker.track(unarchivedFullMetadata)
    _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, unarchivedFullMetadata)
    _ <- waitForArchivedStatus(submittedWorkflow, workflowDefinition)
    workflowEndTime <- extractWorkflowEndTime(unarchivedFullMetadata).asTest

    // Test that these endpoints are still good after archival - but in the window where deletion is unlikely to have happened yet:
    labelsAfterArchival = CentaurCromwellClient.labels(submittedWorkflow)
    _ <- assertLabelsResponseHasNotChanged("after archival but probably not yet deleted", labelsLikelyBeforeArchival, labelsAfterArchival, submittedWorkflow, workflowDefinition)
    _ <- assertMetadataResponseHasNotChanged("after archival but probably not yet deleted", unarchivedNonSubworkflowMetadata, submittedWorkflow, workflowDefinition, expandSubworkflows = false)
    _ <- assertMetadataResponseHasNotChanged("after archival but probably not yet deleted", unarchivedFullMetadata, submittedWorkflow, workflowDefinition, expandSubworkflows = true)
    _ <- assertJobManagerStyleMetadataDidNotChange("after archival but probably not yet deleted", unarchivedJobManagerStyleMetadata, submittedWorkflow, workflowDefinition)
    _ <- assertOutputsResponseHasNotChanged("after archival but probably not yet deleted", unarchivedOutputs, submittedWorkflow, workflowDefinition)

    _ <- waitForArchivedAndPurgedStatus(submittedWorkflow, workflowDefinition, workflowEndTime)

    // Test that these are still good after archival AND deletion:
    labelsAfterArchivalAndDeletion = CentaurCromwellClient.labels(submittedWorkflow)
    _ <- assertLabelsResponseHasNotChanged("after archival and deletion", labelsLikelyBeforeArchival, labelsAfterArchivalAndDeletion, submittedWorkflow, workflowDefinition)
    _ <- assertMetadataResponseHasNotChanged("after archival and deletion", unarchivedNonSubworkflowMetadata, submittedWorkflow, workflowDefinition, expandSubworkflows = false)
    _ <- assertMetadataResponseHasNotChanged("after archival and deletion", unarchivedFullMetadata, submittedWorkflow, workflowDefinition, expandSubworkflows = true)
    _ <- assertJobManagerStyleMetadataDidNotChange("after archival and deletion", unarchivedJobManagerStyleMetadata, submittedWorkflow, workflowDefinition)
    _ <- assertOutputsResponseHasNotChanged("after archival and deletion", unarchivedOutputs, submittedWorkflow, workflowDefinition)

    // Now that we've checked things haven't changed unexpectedly, let's check that we _can_ still change the labels:
    _ <- validateLabelsAdditionUpdateAndSubsequentRetrieval(submittedWorkflow, workflowDefinition)

    // validateLogs method validates data it gets from the `logs` endpoint against logs extracted from the provided
    // metadata (notArchivedMetadata in this case), so there is no need for direct old-new comparison
    _ <- validateLogs(unarchivedNonSubworkflowMetadata, submittedWorkflow, workflowDefinition, validateArchived = Option(true))
  } yield SubmitResponse(submittedWorkflow)

  def runFailingWorkflowAndVerifyMetadata(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = for {
    _ <- checkDescription(workflowDefinition, validityExpectation = None)
    submittedWorkflow <- runFailingWorkflow(workflowDefinition)
    labelsLikelyBeforeArchival = CentaurCromwellClient.labels(submittedWorkflow)
    unarchivedNonSubworkflowMetadata <- fetchAndValidateNonSubworkflowMetadata(submittedWorkflow, workflowDefinition, validateArchived = Option(false))
    unarchivedFullMetadata <- fetchMetadata(submittedWorkflow, expandSubworkflows = true, requestArchivedMetadata = Option(false)).asTest
    unarchivedJobManagerStyleMetadata <- fetchAndValidateJobManagerStyleMetadata(submittedWorkflow, workflowDefinition, prefetchedOriginalNonSubWorkflowMetadata = Option(unarchivedNonSubworkflowMetadata.value), validateArchived = Option(false))
    workflowRoot = unarchivedFullMetadata.asFlat.value.get("workflowRoot").collectFirst { case JsString(r) => r } getOrElse "No Workflow Root"
    unarchivedOutputs <- fetchAndValidateOutputs(submittedWorkflow, workflowDefinition, workflowRoot, validateArchived = Option(false))
    _ = cromwellTracker.track(unarchivedFullMetadata)
    _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, unarchivedFullMetadata)
    _ <- waitForArchivedStatus(submittedWorkflow, workflowDefinition)
    workflowEndTime <- extractWorkflowEndTime(unarchivedFullMetadata).asTest

    _ <- waitForArchivedAndPurgedStatus(submittedWorkflow, workflowDefinition, workflowEndTime)
    // Re-validate the metadata now that carboniting has completed
    labelsAfterArchivalAndDeletion = CentaurCromwellClient.labels(submittedWorkflow)
    _ <- assertLabelsResponseHasNotChanged("after archival and deletion", labelsLikelyBeforeArchival, labelsAfterArchivalAndDeletion, submittedWorkflow, workflowDefinition)
    _ <- assertMetadataResponseHasNotChanged("after archival and deletion", unarchivedNonSubworkflowMetadata, submittedWorkflow, workflowDefinition, expandSubworkflows = false)
    _ <- assertMetadataResponseHasNotChanged("after archival and deletion", unarchivedFullMetadata, submittedWorkflow, workflowDefinition, expandSubworkflows = true)
    _ <- assertJobManagerStyleMetadataDidNotChange("after archival and deletion", unarchivedJobManagerStyleMetadata, submittedWorkflow, workflowDefinition)
    _ <- assertOutputsResponseHasNotChanged("after archival and deletion", unarchivedOutputs, submittedWorkflow, workflowDefinition)

    // Now that we've checked things haven't changed unexpectedly, let's check that we _can_ still change the labels:
    _ <- validateLabelsAdditionUpdateAndSubsequentRetrieval(submittedWorkflow, workflowDefinition)
  } yield SubmitResponse(submittedWorkflow)

  def runWorkflowTwiceExpectingCaching(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    for {
      _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
      firstWF <- runSuccessfulWorkflow(workflowDefinition)
      secondWf <- runSuccessfulWorkflow(workflowDefinition.secondRun)
      _ <- printHashDifferential(firstWF, secondWf)
      metadata <- fetchAndValidateNonSubworkflowMetadata(secondWf, workflowDefinition, Option(firstWF.id.id))
      _ <- fetchAndValidateJobManagerStyleMetadata(secondWf, workflowDefinition, prefetchedOriginalNonSubWorkflowMetadata = None)
      _ = cromwellTracker.track(metadata)
      _ <- validateNoCacheMisses(secondWf, metadata, workflowDefinition)
      _ <- validateDirectoryContentsCounts(workflowDefinition, secondWf, metadata)
    } yield SubmitResponse(secondWf)
  }

  def runWorkflowThriceExpectingCaching(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    for {
      _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
      firstWf <- runSuccessfulWorkflow(workflowDefinition)
      secondWf <- runSuccessfulWorkflow(workflowDefinition.secondRun)
      metadataTwo <- fetchAndValidateNonSubworkflowMetadata(secondWf, workflowDefinition, Option(firstWf.id.id))
      _ = cromwellTracker.track(metadataTwo)
      _ <- validateNoCacheHits(secondWf, metadataTwo, workflowDefinition)
      thirdWf <- runSuccessfulWorkflow(workflowDefinition.thirdRun)
      _ <- printHashDifferential(secondWf, thirdWf)
      metadataThree <- fetchAndValidateNonSubworkflowMetadata(thirdWf, workflowDefinition, Option(secondWf.id.id))
      _ <- validateNoCacheMisses(thirdWf, metadataThree, workflowDefinition)
      _ <- validateDirectoryContentsCounts(workflowDefinition, thirdWf, metadataThree)
    } yield SubmitResponse(thirdWf)
  }

  def runWorkflowTwiceExpectingNoCaching(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    for {
      _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
      _ <- runSuccessfulWorkflow(workflowDefinition) // Build caches
      testWf <- runSuccessfulWorkflow(workflowDefinition.secondRun)
      metadata <- fetchAndValidateNonSubworkflowMetadata(testWf, workflowDefinition)
      _ <- fetchAndValidateJobManagerStyleMetadata(testWf, workflowDefinition, prefetchedOriginalNonSubWorkflowMetadata = None)
      _ = cromwellTracker.track(metadata)
      _ <- validateNoCacheHits(testWf, metadata, workflowDefinition)
      _ <- validateDirectoryContentsCounts(workflowDefinition, testWf, metadata)
    } yield SubmitResponse(testWf)
  }

  def runFailingWorkflowTwiceExpectingNoCaching(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    for {
      _ <- checkDescription(workflowDefinition, validityExpectation = None)
      _ <- runFailingWorkflow(workflowDefinition) // Build caches
      testWf <- runFailingWorkflow(workflowDefinition)
      metadata <- fetchAndValidateNonSubworkflowMetadata(testWf, workflowDefinition)
      _ <- fetchAndValidateJobManagerStyleMetadata(testWf, workflowDefinition, prefetchedOriginalNonSubWorkflowMetadata = None)
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
          _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
          submittedWorkflow <- submitWorkflow(workflowDefinition)
          jobId <- pollUntilCallIsRunning(workflowDefinition, submittedWorkflow, callMarker.callKey)
          _ = CromwellManager.stopCromwell(s"Scheduled restart from ${workflowDefinition.testName}")
          _ = CromwellManager.startCromwell(postRestart)
          _ <- expectSomeProgress(submittedWorkflow, workflowDefinition, Set(Running, finalStatus), workflowProgressTimeout)
          _ <- pollUntilStatus(submittedWorkflow, workflowDefinition, finalStatus)
          metadata <- fetchAndValidateNonSubworkflowMetadata(submittedWorkflow, workflowDefinition)
          _ <- fetchAndValidateJobManagerStyleMetadata(submittedWorkflow, workflowDefinition, prefetchedOriginalNonSubWorkflowMetadata = None)
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
    _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
    submittedWorkflow <- submitWorkflow(workflowDefinition)
    _ <- abortWorkflow(submittedWorkflow)
    _ <- expectSomeProgress(submittedWorkflow, workflowDefinition, Set(Running, Aborting, Aborted), workflowProgressTimeout)
    _ <- pollUntilStatus(submittedWorkflow, workflowDefinition, Aborted)
    metadata <- fetchAndValidateNonSubworkflowMetadata(submittedWorkflow, workflowDefinition)
    _ <- fetchAndValidateJobManagerStyleMetadata(submittedWorkflow, workflowDefinition, prefetchedOriginalNonSubWorkflowMetadata = None)
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
      _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
      submittedWorkflow <- submitWorkflow(workflowDefinition)
      jobId <- pollUntilCallIsRunning(workflowDefinition, submittedWorkflow, callMarker.callKey)
      // The Cromwell call status could be running but the backend job might not have started yet, give it some time
      _ <- waitFor(30.seconds)
      _ <- abortWorkflow(submittedWorkflow)
      _ = if(restart) withRestart()
      _ <- expectSomeProgress(submittedWorkflow, workflowDefinition, Set(Running, Aborting, Aborted), workflowProgressTimeout)
      _ <- pollUntilStatus(submittedWorkflow, workflowDefinition, Aborted)
      _ <- validatePAPIAborted(workflowDefinition, submittedWorkflow, jobId)
      // Wait a little to make sure that if the abort didn't work and calls start running we see them in the metadata
      _ <- waitFor(30.seconds)
      metadata <- fetchAndValidateNonSubworkflowMetadata(submittedWorkflow, workflowDefinition)
      _ <- fetchAndValidateJobManagerStyleMetadata(submittedWorkflow, workflowDefinition, prefetchedOriginalNonSubWorkflowMetadata = None)
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
      _ <- checkDescription(workflow, validityExpectation = None)
      actualSubmitResponse <- Operations.submitInvalidWorkflow(workflow)
      _ <- validateSubmitFailure(workflow, expectedSubmitResponse, actualSubmitResponse)
    } yield actualSubmitResponse
  }

  def papiUpgrade(workflowDefinition: Workflow, callMarker: CallMarker)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    CentaurConfig.runMode match {
      case ManagedCromwellServer(_, postRestart, withRestart) if withRestart =>
        for {
          _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
          first <- submitWorkflow(workflowDefinition)
          _ <- pollUntilCallIsRunning(workflowDefinition, first, callMarker.callKey)
          _ = CromwellManager.stopCromwell(s"Scheduled restart from ${workflowDefinition.testName}")
          _ = CromwellManager.startCromwell(postRestart)
          _ <- expectSomeProgress(first, workflowDefinition, Set(Running, Succeeded), workflowProgressTimeout)
          _ <- pollUntilStatus(first, workflowDefinition, Succeeded)
          _ <- checkDescription(workflowDefinition.secondRun, validityExpectation = Option(true))
          second <- runSuccessfulWorkflow(workflowDefinition.secondRun) // Same WDL and config but a "backend" runtime option targeting PAPI v2.
          _ <- printHashDifferential(first, second)
          metadata <- fetchAndValidateNonSubworkflowMetadata(second, workflowDefinition, Option(first.id.id))
          _ <- fetchAndValidateJobManagerStyleMetadata(second, workflowDefinition, prefetchedOriginalNonSubWorkflowMetadata = None)
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
