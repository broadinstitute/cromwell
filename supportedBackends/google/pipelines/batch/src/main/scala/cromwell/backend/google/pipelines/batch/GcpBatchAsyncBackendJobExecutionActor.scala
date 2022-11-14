package cromwell.backend.google.pipelines.batch

import cromwell.backend.standard.{StandardAsyncExecutionActor, StandardAsyncExecutionActorParams}
import cromwell.core.retry.SimpleExponentialBackoff
import scala.concurrent.duration._
import cromwell.backend._


class GcpBatchAsyncBackendJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams) extends BackendJobLifecycleActor with StandardAsyncExecutionActor {
  /** The type of the run info when a job is started. */
  override type StandardAsyncRunInfo = this.type
  /** The type of the run status returned during each poll. */
  override type StandardAsyncRunState = this.type

  /** Should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be carried around
    * in the state type.
    */
  override def statusEquivalentTo(thiz: GcpBatchAsyncBackendJobExecutionActor.this.type)(that: GcpBatchAsyncBackendJobExecutionActor.this.type): Boolean = true

  /**
    * Returns true when a job is complete, either successfully or unsuccessfully.
    *
    * @param runStatus The run status.
    * @return True if the job has completed.
    */
  override def isTerminal(runStatus: GcpBatchAsyncBackendJobExecutionActor.this.type): Boolean = true

  override def dockerImageUsed: Option[String] = Option("test")


  //override def pollBackOff: SimpleExponentialBackoff = ???

  override lazy val pollBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(1.second, 5.minutes, 1.1)

  //override def executeOrRecoverBackOff: SimpleExponentialBackoff = ???

  override lazy val executeOrRecoverBackOff: SimpleExponentialBackoff = SimpleExponentialBackoff(
    initialInterval = 3.second, maxInterval = 20.second, multiplier = 1.1)
}
