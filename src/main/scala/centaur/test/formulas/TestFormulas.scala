package centaur.test.formulas

import cats.syntax.flatMap._
import cats.syntax.functor._
import centaur.test.Operations._
import centaur.test.Test.testMonad
import centaur.test.markers.CallMarker
import centaur.test.submit.SubmitResponse
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
      _ <- pollUntilStatus(s, status)
    } yield s
  }

  private def runSuccessfulWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = runWorkflowUntilTerminalStatus(workflow, Succeeded)
  private def runFailingWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = runWorkflowUntilTerminalStatus(workflow, Failed)

  def runSuccessfulWorkflowAndVerifyMetadata(workflowDefinition: Workflow): Test[Unit] = for {
    w <- runSuccessfulWorkflow(workflowDefinition)
    _ <- validateMetadata(w, workflowDefinition)
    _ <- validateDirectoryContentsCounts(workflowDefinition, w)
  } yield ()

  def runFailingWorkflowAndVerifyMetadata(workflowDefinition: Workflow): Test[Unit] = for {
    w <- runFailingWorkflow(workflowDefinition)
    _ <- validateMetadata(w, workflowDefinition)
    _ <- validateDirectoryContentsCounts(workflowDefinition, w)
  } yield ()

  def runWorkflowTwiceExpectingCaching(workflowDefinition: Workflow): Test[Unit] = {
    for {
      firstWF <- runSuccessfulWorkflow(workflowDefinition)
      secondWf <- runSuccessfulWorkflow(workflowDefinition)
      metadata <- validateMetadata(secondWf, workflowDefinition, Option(firstWF.id.id))
      _ <- validateNoCacheMisses(metadata, workflowDefinition.testName)
      _ <- validateDirectoryContentsCounts(workflowDefinition, secondWf)
    } yield ()
  }

  def runWorkflowTwiceExpectingNoCaching(workflowDefinition: Workflow): Test[Unit] = {
    for {
      _ <- runSuccessfulWorkflow(workflowDefinition) // Build caches
      testWf <- runSuccessfulWorkflow(workflowDefinition)
      metadata <- validateMetadata(testWf, workflowDefinition)
      _ <- validateNoCacheHits(metadata, workflowDefinition.testName)
      _ <- validateDirectoryContentsCounts(workflowDefinition, testWf)
    } yield ()
  }

  def runFailingWorkflowTwiceExpectingNoCaching(workflowDefinition: Workflow): Test[Unit] = {
    for {
      _ <- runFailingWorkflow(workflowDefinition) // Build caches
      testWf <- runFailingWorkflow(workflowDefinition)
      metadata <- validateMetadata(testWf, workflowDefinition)
      _ <- validateNoCacheHits(metadata, workflowDefinition.testName)
      _ <- validateDirectoryContentsCounts(workflowDefinition, testWf)
    } yield ()
  }
  
  private def cromwellRestart(workflowDefinition: Workflow, callMarker: CallMarker, testRecover: Boolean): Test[Unit] = CentaurConfig.runMode match {
    case ManagedCromwellServer(_, postRestart, withRestart) if withRestart =>
      for {
        w <- submitWorkflow(workflowDefinition)
        jobId <- pollUntilCallIsRunning(w, callMarker.callKey)
        _ = CromwellManager.stopCromwell()
        _ = CromwellManager.startCromwell(postRestart)
        _ <- pollUntilStatus(w, Succeeded)
        _ <- validateMetadata(w, workflowDefinition)
        _ <- if(testRecover) validateRecovered(w, callMarker.callKey, jobId) else Test.successful(())
        _ <- validateDirectoryContentsCounts(workflowDefinition, w)
      } yield ()
    case _ => runSuccessfulWorkflowAndVerifyMetadata(workflowDefinition)
  }

  def instantAbort(workflowDefinition: Workflow): Test[Unit] = for {
    w <- submitWorkflow(workflowDefinition)
    _ <- abortWorkflow(w)
    _ <- pollUntilStatus(w, Aborted)
    _ <- validateMetadata(w, workflowDefinition)
    _ <- validateDirectoryContentsCounts(workflowDefinition, w)
  } yield ()

  def scheduledAbort(workflowDefinition: Workflow, callMarker: CallMarker, restart: Boolean): Test[Unit] = {
    def withRestart() = CentaurConfig.runMode match {
      case ManagedCromwellServer(_, postRestart, withRestart) if withRestart =>
        CromwellManager.stopCromwell()
        CromwellManager.startCromwell(postRestart)
      case _ =>
    }
    
    for {
      w <- submitWorkflow(workflowDefinition)
      jobId <- pollUntilCallIsRunning(w, callMarker.callKey)
      _ <- abortWorkflow(w)
      _ = if(restart) withRestart()
      _ <- pollUntilStatus(w, Aborted)
      _ <- validatePAPIAborted(jobId, w)
      // Wait a little to make sure that if the abort didn't work and calls start running we see them in the metadata
      _ <- waitFor(30.seconds)
      _ <- validateMetadata(w, workflowDefinition)
      _ <- validateDirectoryContentsCounts(workflowDefinition, w)
    } yield ()
  }

  def cromwellRestartWithRecover(workflowDefinition: Workflow, callMarker: CallMarker): Test[Unit] = {
    cromwellRestart(workflowDefinition, callMarker, testRecover = true)
  }

  def cromwellRestartWithoutRecover(workflowDefinition: Workflow, callMarker: CallMarker): Test[Unit] = {
    cromwellRestart(workflowDefinition, callMarker, testRecover = false)
  }

  def submitInvalidWorkflow(workflow: Workflow, expectedSubmitResponse: SubmitResponse): Test[Unit] = {
    for {
      actualSubmitResponse <- Operations.submitInvalidWorkflow(workflow)
      _ <- validateSubmitFailure(workflow, expectedSubmitResponse, actualSubmitResponse)
    } yield ()
  }
}
