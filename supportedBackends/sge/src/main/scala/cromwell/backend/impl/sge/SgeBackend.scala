package cromwell.backend.impl.sge

import cromwell.backend.WorkflowBackendActor.{AbortResponse, ExecutionResponse}
import cromwell.backend._

import scala.concurrent.{ExecutionContext, Future}

class SgeBackend extends WorkflowBackendActor {
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
