package centaur

import java.nio.file.Path

import cats.implicits._
import Operations._
import Test.testMonad

/**
  * A collection of test formulas which can be used, building upon operations by chaining them together via a
  * for comprehension. These assembled formulas can then be run by a client
  */
object CacheFormulas {

  def runCachingWorkflow(request: WorkflowRequest) = {
    for {
      testWF <- TestFormulas.runWorkflowUntilTerminalStatus(request, Succeeded)
      _ <- verifyInputsOutputs(testWF, request)
      _ <- verifyCaching(testWF, request)
    } yield ()
  }

  def runSequentialCachingWorkflow(request: WorkflowRequest, chainedRequest: WorkflowRequest) = {
    for {
      _ <- TestFormulas.runWorkflowUntilTerminalStatus(request, Succeeded)
      _ <- runCachingWorkflow(chainedRequest)
    } yield ()
  }

  def runCachingTurnedOffWorkflow(request: WorkflowRequest) = {
    for {
      testWF <- TestFormulas.runWorkflowUntilTerminalStatus(request, Succeeded)
      _ <- verifyCachingOff(testWF, request)
    } yield ()
  }

}
