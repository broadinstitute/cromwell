package cromwell.engine.callexecution

import com.google.api.client.util.ExponentialBackOff
import cromwell.engine.backend._
import cromwell.engine.callexecution.CallExecutionActor._
import cromwell.logging.WorkflowLogger

import scala.concurrent.ExecutionContext
import scala.language.postfixOps

/**
  * Actor to manage the execution of a single backend call.
  * */
class BackendCallExecutionActor(jobDescriptor: BackendCallJobDescriptor) extends CallExecutionActor {
  override val logger = WorkflowLogger(
    this.getClass.getSimpleName,
    jobDescriptor.workflowDescriptor,
    akkaLogger = Option(akkaLogger),
    callTag = Option(jobDescriptor.key.tag)
  )

  override val call = jobDescriptor.key.scope
  override def poll(handle: ExecutionHandle) = jobDescriptor.poll(handle)
  override def execute(mode: ExecutionMode)(implicit ec: ExecutionContext) = mode match {
    case Execute => jobDescriptor.execute
    case Resume(jobKey) => jobDescriptor.resume(jobKey)
    case UseCachedCall(cachedBackendCall) => jobDescriptor.useCachedCall(cachedBackendCall)
  }

  override val backoff: ExponentialBackOff = jobDescriptor.backend.pollBackoff
}
