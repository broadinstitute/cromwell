package cromwell.backend.impl.sge

import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.BackendLifecycleActor.JobAbortResponse
import cromwell.backend._
import wdl4s.Call

import scala.concurrent.{ExecutionContext, Future}

class SgeJobExecutionActor extends BackendJobExecutionActor {
  override protected implicit def ec: ExecutionContext = ???

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
  override def abortJob: Future[JobAbortResponse] = ???
}
