package cromwell.backend.impl.jes

import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.BackendLifecycleActor.JobAbortResponse
import cromwell.backend._
import wdl4s.Call

import scala.concurrent.Future

class JesBackend extends BackendJobExecutionActor {
  /**
    * Restart or resume a previously-started job.
    */
  override def recover: Future[BackendJobExecutionResponse] = ???

  /**
    * Execute a new job.
    */
  override def execute: Future[BackendJobExecutionResponse] = ???

  override protected def configurationDescriptor: BackendConfigurationDescriptor = ???

  override protected def jobDescriptor: BackendJobDescriptor = ???

  /**
    * Abort a running job.
    */
  override def abortJob: Future[JobAbortResponse] = ???
}
