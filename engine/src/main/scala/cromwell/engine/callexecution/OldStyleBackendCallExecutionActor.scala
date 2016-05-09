package cromwell.engine.callexecution

import com.google.api.client.util.ExponentialBackOff
import cromwell.engine.backend._
import cromwell.engine.callexecution.OldStyleCallExecutionActor._
import cromwell.logging.WorkflowLogger

import scala.concurrent.ExecutionContext
import scala.language.postfixOps

/**
  * Actor to manage the execution of a single backend call.
  * */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class OldStyleBackendCallExecutionActor(jobDescriptor: OldStyleBackendCallJobDescriptor) extends OldStyleCallExecutionActor {
  override val logger = WorkflowLogger(
    this.getClass.getSimpleName,
    jobDescriptor.workflowDescriptor,
    akkaLogger = Option(akkaLogger),
    callTag = Option(jobDescriptor.key.tag)
  )

  override val call = jobDescriptor.call
  override def poll(handle: OldStyleExecutionHandle) = jobDescriptor.poll(handle)
  override def execute(mode: ExecutionMode)(implicit ec: ExecutionContext) = mode match {
    case Execute => jobDescriptor.execute
    case Resume(executionInfos) => jobDescriptor.resume(executionInfos)
    case UseCachedCall(cachedBackendCall) => jobDescriptor.useCachedCall(cachedBackendCall)
  }

  override val backoff: ExponentialBackOff = jobDescriptor.backend.pollBackoff
}
