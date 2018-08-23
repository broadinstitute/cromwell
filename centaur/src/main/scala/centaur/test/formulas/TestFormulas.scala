package centaur.test.formulas

import cats.syntax.flatMap._
import cats.syntax.functor._
import centaur.test.Operations._
import centaur.test.Test.testMonad
import centaur.test.markers.CallMarker
import centaur.test.submit.{SubmitHttpResponse, SubmitResponse}
import centaur.test.workflow.Workflow
import centaur.test.{Operations, Test}
import centaur.{CentaurConfig, CromwellManager, ManagedCromwellServer}
import cromwell.api.model.{Aborted, Failed, SubmittedWorkflow, Succeeded, TerminalStatus}

import scala.concurrent.duration._

/**
  * A collection of test formulas which can be used, building upon operations by chaining them together via a
  * for comprehension. These assembled formulas can then be run by a client
  */
object TestFormulas {
  private def runWorkflowUntilTerminalStatus(workflow: Workflow, status: TerminalStatus): Test[SubmittedWorkflow] = {
    for {
      s <- submitWorkflow(workflow)
      _ <- pollUntilStatus(s, workflow, status)
    } yield s
  }

  private def runSuccessfulWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = runWorkflowUntilTerminalStatus(workflow, Succeeded)
  private def runFailingWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = runWorkflowUntilTerminalStatus(workflow, Failed)

  def runSuccessfulWorkflowAndVerifyMetadata(workflowDefinition: Workflow): Test[SubmitResponse] = for {
    submittedWorkflow <- runSuccessfulWorkflow(workflowDefinition)
    metadata <- validateMetadata(submittedWorkflow, workflowDefinition)
    _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, metadata)
  } yield SubmitResponse(submittedWorkflow)

  def runFailingWorkflowAndVerifyMetadata(workflowDefinition: Workflow): Test[SubmitResponse] = for {
    submittedWorkflow <- runFailingWorkflow(workflowDefinition)
    metadata <- validateMetadata(submittedWorkflow, workflowDefinition)
    _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, metadata)
  } yield SubmitResponse(submittedWorkflow)

  def runWorkflowTwiceExpectingCaching(workflowDefinition: Workflow): Test[SubmitResponse] = {
    for {
      firstWF <- runSuccessfulWorkflow(workflowDefinition)
      secondWf <- runSuccessfulWorkflow(workflowDefinition)
      _ <- printHashDifferential(firstWF, secondWf)
      metadata <- validateMetadata(secondWf, workflowDefinition, Option(firstWF.id.id))
      _ <- validateNoCacheMisses(secondWf, metadata, workflowDefinition)
      _ <- validateDirectoryContentsCounts(workflowDefinition, secondWf, metadata)
    } yield SubmitResponse(secondWf)
  }

  def runWorkflowTwiceExpectingNoCaching(workflowDefinition: Workflow): Test[SubmitResponse] = {
    for {
      _ <- runSuccessfulWorkflow(workflowDefinition) // Build caches
      testWf <- runSuccessfulWorkflow(workflowDefinition)
      metadata <- validateMetadata(testWf, workflowDefinition)
      _ <- validateNoCacheHits(testWf, metadata, workflowDefinition)
      _ <- validateDirectoryContentsCounts(workflowDefinition, testWf, metadata)
    } yield SubmitResponse(testWf)
  }

  def runFailingWorkflowTwiceExpectingNoCaching(workflowDefinition: Workflow): Test[SubmitResponse] = {
    for {
      _ <- runFailingWorkflow(workflowDefinition) // Build caches
      testWf <- runFailingWorkflow(workflowDefinition)
      metadata <- validateMetadata(testWf, workflowDefinition)
      _ <- validateNoCacheHits(testWf, metadata, workflowDefinition)
      _ <- validateDirectoryContentsCounts(workflowDefinition, testWf, metadata)
    } yield SubmitResponse(testWf)
  }
  
  private def cromwellRestart(workflowDefinition: Workflow,
                              callMarker: CallMarker,
                              testRecover: Boolean,
                              finalStatus: TerminalStatus): Test[SubmitResponse] = {
    CentaurConfig.runMode match {
      case ManagedCromwellServer(_, postRestart, withRestart) if withRestart =>
        for {
          submittedWorkflow <- submitWorkflow(workflowDefinition)
          jobId <- pollUntilCallIsRunning(workflowDefinition, submittedWorkflow, callMarker.callKey)
          _ = CromwellManager.stopCromwell(s"Scheduled restart from ${workflowDefinition.testName}")
          _ = CromwellManager.startCromwell(postRestart)
          _ <- pollUntilStatus(submittedWorkflow, workflowDefinition, finalStatus)
          metadata <- validateMetadata(submittedWorkflow, workflowDefinition)
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

  def instantAbort(workflowDefinition: Workflow): Test[SubmitResponse] = for {
    submittedWorkflow <- submitWorkflow(workflowDefinition)
    _ <- abortWorkflow(submittedWorkflow)
    _ <- pollUntilStatus(submittedWorkflow, workflowDefinition, Aborted)
    metadata <- validateMetadata(submittedWorkflow, workflowDefinition)
    _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, metadata)
  } yield SubmitResponse(submittedWorkflow)

  def scheduledAbort(workflowDefinition: Workflow, callMarker: CallMarker, restart: Boolean): Test[SubmitResponse] = {
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
      _ <- pollUntilStatus(submittedWorkflow, workflowDefinition, Aborted)
      _ <- validatePAPIAborted(workflowDefinition, submittedWorkflow, jobId)
      // Wait a little to make sure that if the abort didn't work and calls start running we see them in the metadata
      _ <- waitFor(30.seconds)
      metadata <- validateMetadata(submittedWorkflow, workflowDefinition)
      _ <- validateDirectoryContentsCounts(workflowDefinition, submittedWorkflow, metadata)
    } yield SubmitResponse(submittedWorkflow)
  }

  def workflowRestart(workflowDefinition: Workflow,
                             callMarker: CallMarker,
                             recover: Boolean,
                             finalStatus: TerminalStatus): Test[SubmitResponse] = {
    cromwellRestart(workflowDefinition, callMarker, testRecover = recover, finalStatus = finalStatus)
  }

  def submitInvalidWorkflow(workflow: Workflow, expectedSubmitResponse: SubmitHttpResponse): Test[SubmitResponse] = {
    for {
      actualSubmitResponse <- Operations.submitInvalidWorkflow(workflow)
      _ <- validateSubmitFailure(workflow, expectedSubmitResponse, actualSubmitResponse)
    } yield actualSubmitResponse
  }
}
