package centaur.test.formulas

import java.time.OffsetDateTime

import cats.syntax.flatMap._
import cats.syntax.functor._
import centaur.test.Operations._
import centaur.test.Test._
import centaur.test.markers.CallMarker
import centaur.test.submit.{SubmitHttpResponse, SubmitResponse}
import centaur.test.workflow.Workflow
import centaur.test.{Operations, Test}
import centaur.{CentaurConfig, CromwellManager, CromwellTracker, ManagedCromwellServer}
import com.typesafe.scalalogging.StrictLogging
import cromwell.api.model.{Aborted, Aborting, Failed, Running, SubmittedWorkflow, Succeeded, TerminalStatus, WorkflowMetadata}

import scala.concurrent.duration._
import centaur.test.metadata.WorkflowFlatMetadata._
import spray.json.JsString
/**
  * A collection of test formulas which can be used, building upon operations by chaining them together via a
  * for comprehension. These assembled formulas can then be run by a client
  */
object TestFormulas extends StrictLogging {
  logger.info(
    s"""Running with CentaurConfig:
       |  - expectCarbonite: '${CentaurConfig.expectCarbonite}'
       |  - workflowProgressTimeout: '${CentaurConfig.workflowProgressTimeout}'
       |  - sendReceiveTimeout: '${CentaurConfig.sendReceiveTimeout}'
       |  - maxWorkflowLength: '${CentaurConfig.maxWorkflowLength}'
       |  - metadataConsistencyTimeout: '${CentaurConfig.metadataConsistencyTimeout}'
       |  - standardTestCasePath: '${CentaurConfig.standardTestCasePath}'
       |  - optionalTestPath: '${CentaurConfig.optionalTestPath}'
       |""".stripMargin.trim
  )

  private def runWorkflowUntilTerminalStatus(workflow: Workflow, status: TerminalStatus): Test[SubmittedWorkflow] = {
    for {
      _ <- checkVersion()
      s <- submitWorkflow(workflow)
      _ <- expectSomeProgress(
        workflow = s,
        testDefinition = workflow,
        expectedStatuses = Set(Running, status),
        timeout = CentaurConfig.workflowProgressTimeout,
      )
      _ <- pollUntilStatus(s, workflow, status)
    } yield s
  }

  private def runSuccessfulWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = runWorkflowUntilTerminalStatus(workflow, Succeeded)
  private def runFailingWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = runWorkflowUntilTerminalStatus(workflow, Failed)

  def runSuccessfulWorkflowAndVerifyTimeAndOutputs(workflowDefinition: Workflow): Test[SubmitResponse] = for {
    _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
    timeAllowance <- checkTimingRequirement(workflowDefinition.maximumAllowedTime)
    beforeTimestamp = OffsetDateTime.now().toInstant.getEpochSecond
    submittedWorkflow <- runSuccessfulWorkflow(workflowDefinition)
    afterTimestamp = OffsetDateTime.now().toInstant.getEpochSecond
    _ <- fetchAndValidateOutputs(submittedWorkflow, workflowDefinition, "ROOT NOT SUPPORTED IN TIMING/OUTPUT ONLY TESTS")
    _ <- checkFastEnough(beforeTimestamp, afterTimestamp, timeAllowance)
  } yield SubmitResponse(submittedWorkflow)

  def runSuccessfulWorkflowAndVerifyMetadata(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = for {
    _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
    _ <- timingVerificationNotSupported(workflowDefinition.maximumAllowedTime)
    submittedWorkflow <- runSuccessfulWorkflow(workflowDefinition)
    metadata <- fetchAndValidateNonSubworkflowMetadata(submittedWorkflow, workflowDefinition)
    _ <- fetchAndValidateJobManagerStyleMetadata(submittedWorkflow, workflowDefinition, prefetchedOriginalNonSubWorkflowMetadata = Option(metadata.value))
    notArchivedFlatMetadata = metadata.asFlat
    workflowRoot = notArchivedFlatMetadata.value.get("workflowRoot").collectFirst { case JsString(r) => r } getOrElse "No Workflow Root"
    _ <- fetchAndValidateOutputs(submittedWorkflow, workflowDefinition, workflowRoot)
    _ <- fetchAndValidateLabels(submittedWorkflow, workflowDefinition, workflowRoot)
    _ <- validateLogs(metadata, submittedWorkflow, workflowDefinition)
    _ = cromwellTracker.track(metadata)
    _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, metadata)
  } yield SubmitResponse(submittedWorkflow)

  def runFailingWorkflowAndVerifyMetadata(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = for {
    _ <- checkDescription(workflowDefinition, validityExpectation = None)
    _ <- timingVerificationNotSupported(workflowDefinition.maximumAllowedTime)
    submittedWorkflow <- runFailingWorkflow(workflowDefinition)
    metadata <- fetchAndValidateNonSubworkflowMetadata(submittedWorkflow, workflowDefinition)
    _ <- fetchAndValidateJobManagerStyleMetadata(submittedWorkflow, workflowDefinition, Option(metadata.value))
    _ = cromwellTracker.track(metadata)
    _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, metadata)
  } yield SubmitResponse(submittedWorkflow)

  def runWorkflowTwiceExpectingCaching(workflowDefinition: Workflow)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    for {
      _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
      _ <- timingVerificationNotSupported(workflowDefinition.maximumAllowedTime)
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
      _ <- timingVerificationNotSupported(workflowDefinition.maximumAllowedTime)
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
      _ <- timingVerificationNotSupported(workflowDefinition.maximumAllowedTime)
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
      _ <- timingVerificationNotSupported(workflowDefinition.maximumAllowedTime)
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
          _ <- timingVerificationNotSupported(workflowDefinition.maximumAllowedTime)
          submittedWorkflow <- submitWorkflow(workflowDefinition)
          jobId <- pollUntilCallIsRunning(workflowDefinition, submittedWorkflow, callMarker.callKey)
          _ = CromwellManager.stopCromwell(s"Scheduled restart from ${workflowDefinition.testName}")
          _ = CromwellManager.startCromwell(postRestart)
          _ <- expectSomeProgress(
            workflow = submittedWorkflow,
            testDefinition = workflowDefinition,
            expectedStatuses = Set(Running, finalStatus),
            timeout = CentaurConfig.workflowProgressTimeout,
          )
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
    _ <- timingVerificationNotSupported(workflowDefinition.maximumAllowedTime)
    submittedWorkflow <- submitWorkflow(workflowDefinition)
    _ <- abortWorkflow(submittedWorkflow)
    _ <- expectSomeProgress(
      workflow = submittedWorkflow,
      testDefinition = workflowDefinition,
      expectedStatuses = Set(Running, Aborting, Aborted),
      timeout = CentaurConfig.workflowProgressTimeout,
    )
    _ <- pollUntilStatus(submittedWorkflow, workflowDefinition, Aborted)
    metadata <- fetchAndValidateNonSubworkflowMetadata(submittedWorkflow, workflowDefinition)
    _ <- fetchAndValidateJobManagerStyleMetadata(submittedWorkflow, workflowDefinition, prefetchedOriginalNonSubWorkflowMetadata = None)
    _ = cromwellTracker.track(metadata)
    _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, metadata)
  } yield SubmitResponse(submittedWorkflow)

  def scheduledAbort(workflowDefinition: Workflow, callMarker: CallMarker, restart: Boolean)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    def withRestart(): Unit = CentaurConfig.runMode match {
      case ManagedCromwellServer(_, postRestart, withRestart) if withRestart =>
        CromwellManager.stopCromwell(s"Scheduled restart from ${workflowDefinition.testName}")
        CromwellManager.startCromwell(postRestart)
      case _ =>
    }

    for {
      _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
      _ <- timingVerificationNotSupported(workflowDefinition.maximumAllowedTime)
      submittedWorkflow <- submitWorkflow(workflowDefinition)
      jobId <- pollUntilCallIsRunning(workflowDefinition, submittedWorkflow, callMarker.callKey)
      // The Cromwell call status could be running but the backend job might not have started yet, give it some time
      _ <- waitFor(30.seconds)
      _ <- abortWorkflow(submittedWorkflow)
      _ = if(restart) withRestart()
      _ <- expectSomeProgress(
        workflow = submittedWorkflow,
        testDefinition = workflowDefinition,
        expectedStatuses = Set(Running, Aborting, Aborted),
        timeout = CentaurConfig.workflowProgressTimeout,
      )
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
      _ <- timingVerificationNotSupported(workflow.maximumAllowedTime)
      actualSubmitResponse <- Operations.submitInvalidWorkflow(workflow)
      _ <- validateSubmitFailure(workflow, expectedSubmitResponse, actualSubmitResponse)
    } yield actualSubmitResponse
  }

  def papiUpgrade(workflowDefinition: Workflow, callMarker: CallMarker)(implicit cromwellTracker: Option[CromwellTracker]): Test[SubmitResponse] = {
    CentaurConfig.runMode match {
      case ManagedCromwellServer(_, postRestart, withRestart) if withRestart =>
        for {
          _ <- checkDescription(workflowDefinition, validityExpectation = Option(true))
          _ <- timingVerificationNotSupported(workflowDefinition.maximumAllowedTime)
          first <- submitWorkflow(workflowDefinition)
          _ <- pollUntilCallIsRunning(workflowDefinition, first, callMarker.callKey)
          _ = CromwellManager.stopCromwell(s"Scheduled restart from ${workflowDefinition.testName}")
          _ = CromwellManager.startCromwell(postRestart)
          _ <- expectSomeProgress(
            workflow = first,
            testDefinition = workflowDefinition,
            expectedStatuses = Set(Running, Succeeded),
            timeout = CentaurConfig.workflowProgressTimeout,
          )
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
