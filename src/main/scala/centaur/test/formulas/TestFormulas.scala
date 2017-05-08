package centaur.test.formulas

import java.util.UUID

import cats.syntax.eq._
import cats.syntax.functor._
import cats.syntax.flatMap._
import cats._
import centaur.test.Operations._
import centaur.test.Test
import centaur.test.Test.testMonad
import centaur.test.workflow.Workflow
import centaur.test.workflow.Workflow.{WorkflowWithMetadata, WorkflowWithoutMetadata}
import cromwell.api.model.{SubmittedWorkflow, Succeeded, Failed, TerminalStatus}
import spray.json._
import io.circe._, io.circe.parser._


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
      _ <- runSuccessfulWorkflowAndVerifyCacheMetadata(secondAttempt, firstWF.id.id)
    } yield ()
  }

  def runCachingTurnedOffWorkflow(workflow: Workflow) = {
    for {
      testWf <- runWorkflowUntilTerminalStatus(workflow, Succeeded)
      metadata <- retrieveMetadata(testWf)
      _ <- validateNoCacheHits(metadata, workflow.testName)
    } yield ()
  }

  def runFinalDirsWorkflow(wf: Workflow, dirOption: String): Test[Unit] =  {
    import spray.json.DefaultJsonProtocol._

    def dirSize = {
      val options = wf.data.options.get
      val outputDir = parse(options).toOption.flatMap(_.\\(dirOption).head.asString).get
      val dir = new java.io.File(outputDir)
      dir.listFiles.size
    }

    for {
      terminatedWf <- runWorkflowUntilTerminalStatus(wf, Succeeded)
      _ = if (dirSize == 0) throw new RuntimeException("no files in output dir!") else ()
    } yield ()
  }
}
