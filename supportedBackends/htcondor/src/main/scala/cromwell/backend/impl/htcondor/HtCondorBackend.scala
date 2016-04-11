package cromwell.backend.impl.htcondor

import cromwell.backend.WorkflowBackendActor.{ExecutionResponse, AbortResponse}
import cromwell.backend._

import scala.concurrent.{ExecutionContext, Future}

class HtCondorBackend extends WorkflowBackendActor {
  override protected implicit def ec: ExecutionContext = ???

  /**
    * Restart or resume a previously-started job.
    */
  override def recover(jobDescriptor: BackendJobDescriptor): Future[ExecutionResponse] = ???

  /**
    * Execute a new job.
    */
  override def execute(jobDescriptor: BackendJobDescriptor): Future[ExecutionResponse] = ???

  override protected def workflowDescriptor: BackendWorkflowDescriptor = ???

  override protected def configurationDescriptor: BackendConfigurationDescriptor = ???

  /**
    * Abort a running job.
    */
  override def abort(jobKey: BackendJobDescriptorKey): Future[AbortResponse] = ???
}
