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
  override def recover(jobDescriptor: BackendJobDescriptor): Future[BackendJobExecutionResponse] = ???

  /**
    * Execute a new job.
    */
  override def execute(jobDescriptor: BackendJobDescriptor): Future[BackendJobExecutionResponse] = ???

  override protected def workflowDescriptor: BackendWorkflowDescriptor = ???

  override protected def configurationDescriptor: BackendConfigurationDescriptor = ???

  /**
    * The set of calls which this backend will be expected to run
    */
  override protected def calls: Seq[Call] = ???

  /**
    * Abort a running job.
    */
  override def abortJob(jobKey: BackendJobDescriptorKey): Future[JobAbortResponse] = ???
}
