package cromwell.backend.impl.htcondor

import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend._

import scala.concurrent.Future

class HtCondorJobExecutionActor extends BackendJobExecutionActor {
  /**
    * Restart or resume a previously-started job.
    */
  override def recover: Future[BackendJobExecutionResponse] = ???

  /**
    * Execute a new job.
    */
  override def execute: Future[BackendJobExecutionResponse] = ???

  override protected def jobDescriptor: BackendJobDescriptor = ???

  override protected def configurationDescriptor: BackendConfigurationDescriptor = ???

  /**
    * Abort a running job.
    */
  override def abortJob: Unit = ???
}
