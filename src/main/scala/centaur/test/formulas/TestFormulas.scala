package centaur.test.formulas

import java.util.UUID

import cats.implicits._
import centaur._
import centaur.test.Operations._
import centaur.test.Test
import centaur.test.Test.testMonad
import centaur.test.workflow.Workflow
import centaur.test.workflow.Workflow.{WorkflowWithMetadata, WorkflowWithoutMetadata}

/**
  * A collection of test formulas which can be used, building upon operations by chaining them together via a
  * for comprehension. These assembled formulas can then be run by a client
  */
object TestFormulas {
  def runWorkflowUntilTerminalStatus(workflow: Workflow, status: TerminalStatus): Test[SubmittedWorkflow] = {
    for {
      s <- submitWorkflow(workflow)
      _ <- pollUntilStatus(s, status)
    } yield s
  }

  def runSuccessfulWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = runWorkflowUntilTerminalStatus(workflow, Succeeded)
  def runFailingWorkflow(workflow: Workflow): Test[SubmittedWorkflow] = runWorkflowUntilTerminalStatus(workflow, Failed)

  def runSuccessfulWorkflowAndVerifyMetadata(workflow: Workflow): Test[Unit] = runWorkflowAndVerifyMetadata(workflow, runSuccessfulWorkflow)
  def runFailingWorkflowAndVerifyMetadata(workflow: Workflow): Test[Unit] = runWorkflowAndVerifyMetadata(workflow, runFailingWorkflow)
  def runSuccessfulWorkflowAndVerifyCacheMetadata(workflow: Workflow, cacheHit: UUID): Test[Unit] = runWorkflowAndVerifyCacheMetadata(workflow, runSuccessfulWorkflow, cacheHit)

  def runWorkflowAndVerifyCacheMetadata(workflow: Workflow, f: Workflow => Test[SubmittedWorkflow], cacheHitWorkflow: UUID): Test[Unit] = {
    workflow match {
      case _: WorkflowWithoutMetadata => throw new Exception("Scala type system: 3, Jeff: 0")
      case r: WorkflowWithMetadata =>
        for {
          w <- f(r)
          _ <- validateMetadata(w, r.metadata, Option(cacheHitWorkflow))
        } yield ()
    }
  }

  def runWorkflowAndVerifyMetadata(workflow: Workflow, f: Workflow => Test[SubmittedWorkflow]): Test[Unit] = {
    // FIXME: This is horrible, but I just wanted to add this and copy/paste was easier than thinking
    workflow match {
      case _: WorkflowWithoutMetadata => throw new Exception("Scala type system: 3, Jeff: 0")
      case r: WorkflowWithMetadata =>
        for {
          w <- f(r)
          _ <- validateMetadata(w, r.metadata)
        } yield ()
    }
  }

  def runSequentialCachingWorkflows(firstAttempt: Workflow, secondAttempt: Workflow) = {
    for {
      firstWF <- runWorkflowUntilTerminalStatus(firstAttempt, Succeeded)
      _ <- runSuccessfulWorkflowAndVerifyCacheMetadata(secondAttempt, firstWF.id)
    } yield ()
  }

  def runCachingTurnedOffWorkflow(workflow: Workflow) = {
    for {
      testWf <- runWorkflowUntilTerminalStatus(workflow, Succeeded)
      metadata <- retrieveMetadata(testWf)
      _ <- validateNoCacheHits(metadata, workflow.name)
    } yield ()
  }
}
