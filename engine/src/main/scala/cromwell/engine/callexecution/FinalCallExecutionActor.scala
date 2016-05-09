package cromwell.engine.callexecution

import com.google.api.client.util.ExponentialBackOff
import com.google.api.client.util.ExponentialBackOff.Builder
import cromwell.engine.backend.OldStyleExecutionHandle
import cromwell.engine.callexecution.OldStyleCallExecutionActor._
import cromwell.engine.finalcall.OldStyleFinalCall
import cromwell.logging.WorkflowLogger
import cromwell.webservice.WorkflowMetadataResponse

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class FinalCallExecutionActor(override val call: OldStyleFinalCall, workflowMetadataResponse: WorkflowMetadataResponse)
  extends OldStyleCallExecutionActor {

  override val logger = WorkflowLogger(
    this.getClass.getSimpleName,
    call.workflow,
    akkaLogger = Option(akkaLogger),
    callTag = Option(call.unqualifiedName)
  )

  override def poll(handle: OldStyleExecutionHandle) = call.poll(ec, handle)
  override def execute(mode: ExecutionMode)(implicit ec: ExecutionContext) = mode match {
    case UseCachedCall(cachedBackendCall) =>
      Future.failed(new UnsupportedOperationException("Cannot use cached results for a FinalCall"))
    case _ => call.execute(workflowMetadataResponse) map { _ => call.handle }
  }

  /** There is no polling for Final Calls (yet) as they are executed in Cromwell. */
  override val backoff: ExponentialBackOff = new Builder()
    .setInitialIntervalMillis(10.seconds.toMillis.toInt)
    .setMaxElapsedTimeMillis(Int.MaxValue)
    .setMaxIntervalMillis(10.minutes.toMillis.toInt)
    .setMultiplier(1.1D)
    .build()
}
