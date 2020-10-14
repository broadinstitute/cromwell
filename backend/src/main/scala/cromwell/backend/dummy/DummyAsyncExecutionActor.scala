package cromwell.backend.dummy

import cromwell.backend.BackendJobLifecycleActor
import cromwell.backend.async.{ExecutionHandle, SuccessfulExecutionHandle}
import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams}
import cromwell.core.CallOutputs
import cromwell.core.retry.SimpleExponentialBackoff

import scala.concurrent.duration._

class DummyAsyncExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends BackendJobLifecycleActor
  with StandardAsyncExecutionActor  {

  /** The type of the run info when a job is started. */
  override type StandardAsyncRunInfo = String
  /** The type of the run status returned during each poll. */
  override type StandardAsyncRunState = String

  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */
  override def statusEquivalentTo(thiz: String)(that: String): Boolean = thiz == that

  /**
    * Returns true when a job is complete, either successfully or unsuccessfully.
    *
    * @param runStatus The run status.
    * @return True if the job has completed.
    */
  override def isTerminal(runStatus: String): Boolean = runStatus == "DummyDone"

  override def dockerImageUsed: Option[String] = None

  override def pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(initialInterval = 1.second, maxInterval = 300.seconds, multiplier = 1.1)

  override def executeOrRecoverBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(initialInterval = 1.second, maxInterval = 300.seconds, multiplier = 1.1)

  override def execute(): ExecutionHandle = SuccessfulExecutionHandle(
    outputs = CallOutputs(Map.empty),
    returnCode = 0,
    jobDetritusFiles = Map.empty,
    executionEvents = Seq.empty,
    resultsClonedFrom = None
  )
}
