package centaur

import java.nio.file.Path

import cats.implicits._
import Operations._
import Test.testMonad

/**
  * A collection of test formulas which can be used, building upon operations by chaining them together via a
  * for comprehension. These assembled formulas can then be run by a client
  */
object TestFormulas {
  def runWorkflowUntilTerminalStatus(request: WorkflowRequest, status: TerminalStatus): Test[Workflow] = {
    for {
      s <- submitWorkflow(request)
      _ <- pollUntilStatus(s, status)
    } yield s
  }

  def runSuccessfulWorkflow(request: WorkflowRequest) = {
    for {
      s <- runWorkflowUntilTerminalStatus(request, Succeeded)
      _ <- verifyMetadataAndOutputs(s, request)
    } yield ()
  }

  def runFailingWorkflow(request: WorkflowRequest) = runWorkflowUntilTerminalStatus(request, Failed)

  def runSubmissionFailureWorkflow(request: WorkflowRequest) = submitWorkflowExpectingRejection(request)
}
