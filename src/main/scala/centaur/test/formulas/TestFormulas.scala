package centaur.test.formulas

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

  def runSuccessfulWorkflowAndVerifyMetadata(workflow: Workflow): Test[Unit] = {
    // FIXME: This is horrible, but I want to get the issue closed, will come back and make this more type safe
    workflow match {
      case _: WorkflowWithoutMetadata => throw new Exception("Scala type system: 1, Jeff: 0")
      case r: WorkflowWithMetadata =>
        for {
          w <- runSuccessfulWorkflow(r)
          m <- retrieveMetadata(w)
          _ <- validateMetadata(m, r.metadata, w.id)
        } yield ()
    }
  }

  def runFailingWorkflowAndVerifyMetadata(workflow: Workflow): Test[Unit] = {
    // FIXME: This is horrible, but I just wanted to add this and copy/paste was easier than thinking
    workflow match {
      case _: WorkflowWithoutMetadata => throw new Exception("Scala type system: 2, Jeff: 0")
      case r: WorkflowWithMetadata =>
        for {
          w <- runFailingWorkflow(r)
          m <- retrieveMetadata(w)
          _ <- validateMetadata(m, r.metadata, w.id)
        } yield ()
    }
  }

  def runSequentialCachingWorkflow(workflow: Workflow, secondWorkflow: Workflow) = {
    for {
      _ <- runWorkflowUntilTerminalStatus(workflow, Succeeded)
      _ <- runSuccessfulWorkflowAndVerifyMetadata(secondWorkflow)
    } yield ()
  }

  def runCachingTurnedOffWorkflow(workflow: Workflow) = {
    for {
      testWf <- runWorkflowUntilTerminalStatus(workflow, Succeeded)
      metadata <- retrieveMetadata(testWf)
      _ <- validateCachingWasOff(metadata, workflow.name)
    } yield ()
  }
}
