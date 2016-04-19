package cromwell.backend.impl.sge

import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.BackendLifecycleActor.JobAbortResponse
import cromwell.backend._
import wdl4s.Call

import scala.concurrent.{ExecutionContext, Future}

class SgeBackend extends BackendJobExecutionActor {
  override protected implicit def ec: ExecutionContext = ???

  /**
    * Restart or resume a previously-started job.
    */
  override def recover(jobDescriptor: BackendJobDescriptor): Future[BackendJobExecutionResponse] = ???

  /**
    * Execute a new job.
    */
  override def execute(jobDescriptor: BackendJobDescriptor): Future[BackendJobExecutionResponse] = ???

  override protected def jobDescriptor: BackendJobDescriptor = ???

  override protected def configurationDescriptor: BackendConfigurationDescriptor = ???

  /**
    * Abort a running job.
    */
  override def abortJob(jobKey: BackendJobDescriptorKey): Future[JobAbortResponse] = ???
}
