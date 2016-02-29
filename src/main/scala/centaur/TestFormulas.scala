package centaur

import centaur.Operation.Test

/**
  * A collection of test formulas which can be used, building upon Operations by chaining them together. The cleanest
  * method of chaining would be a for comprehension
  *
  *  FIXME: The need for the .perform is pretty hokey at the moment but it's a bit of a hack to get things in
  *  the direction I wanted in a faster fashion
  */
object TestFormulas {
  def runWorkflowUntilTerminalStatus(request: WorkflowRequest, status: TerminalStatus): Test = {
    for {
      s <- SubmitWorkflow(request).perform
      _ <- PollUntilTerminal(s, status).perform
    } yield s
  }

  def runSuccessfulWorkflow(request: WorkflowRequest): Test = runWorkflowUntilTerminalStatus(request, Succeeded)
  def runFailingWorkflow(request: WorkflowRequest): Test = runWorkflowUntilTerminalStatus(request, Failed)
}
