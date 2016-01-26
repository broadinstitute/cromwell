package cromwell.engine.callexecution

import com.google.api.client.util.ExponentialBackOff
import com.google.api.client.util.ExponentialBackOff.Builder
import cromwell.engine.backend.ExecutionHandle
import cromwell.engine.callexecution.CallExecutionActor._
import cromwell.engine.finalcall.FinalCall
import cromwell.logging.WorkflowLogger

import scala.concurrent.Future
import scala.concurrent.duration._

case class FinalCallExecutionActor(override val call: FinalCall) extends CallExecutionActor {

  override val logger = WorkflowLogger(
    this.getClass.getSimpleName,
    call.workflow,
    akkaLogger = Option(akkaLogger),
    callTag = Option(call.unqualifiedName)
  )

  override def poll(handle: ExecutionHandle) = call.poll(ec, handle)
  override def execute(mode: ExecutionMode) = mode match {
    case UseCachedCall(cachedBackendCall) => Future.failed(new UnsupportedOperationException("Cannot use cached results for a FinalCall"))
    case _ => call.execute
  }

  /** There is no polling for Final Calls (yet) as they are executed in Cromwell. */
  override val backoff: ExponentialBackOff = new Builder()
    .setInitialIntervalMillis(10.seconds.toMillis.toInt)
    .setMaxElapsedTimeMillis(Int.MaxValue)
    .setMaxIntervalMillis(10.minutes.toMillis.toInt)
    .setMultiplier(1.1D)
    .build()
}
