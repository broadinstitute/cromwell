package centaur
import cats.implicits._
import Operations._
import Test.testMonad

/**
  * A collection of test formulas which can be used, building upon operations by chaining them together via a
  * for comprehension. These assembled formulas can then be run by a client
  */
object TestFormulas {
  def runWorkflowUntilTerminalStatus(request: WorkflowRequest, status: TerminalStatus): Test[Unit] = {
    for {
      s <- submitWorkflow(request)
      _ <- pollUntilStatus(s, status)
    } yield ()
  }

  def runSuccessfulWorkflow(request: WorkflowRequest) = runWorkflowUntilTerminalStatus(request, Succeeded)
  def runFailingWorkflow(request: WorkflowRequest) = runWorkflowUntilTerminalStatus(request, Failed)
}
