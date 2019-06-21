package cromwell.engine.workflow.workflowstore

import java.time.{OffsetDateTime, Duration => JDuration}
import java.util.concurrent.TimeUnit

import akka.actor.{ActorRef, CoordinatedShutdown, Props}
import cats.data.{NonEmptyList, NonEmptyVector}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowId
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.engine.CromwellTerminator
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.WorkflowStoreWriteHeartbeatCommand
import cromwell.services.EnhancedBatchActor
import mouse.all._

import scala.concurrent.Future
import scala.concurrent.duration.FiniteDuration
import scala.util.{Failure, Success, Try}

case class WorkflowStoreHeartbeatWriteActor(workflowStoreAccess: WorkflowStoreAccess,
                                            workflowHeartbeatConfig: WorkflowHeartbeatConfig,
                                            terminator: CromwellTerminator,
                                            override val serviceRegistryActor: ActorRef)

  extends EnhancedBatchActor[(WorkflowId, OffsetDateTime)](
    flushRate = workflowHeartbeatConfig.heartbeatInterval,
    batchSize = workflowHeartbeatConfig.writeBatchSize) {

  override val threshold = workflowHeartbeatConfig.writeThreshold

  private val failureShutdownDuration = workflowHeartbeatConfig.failureShutdownDuration

  //noinspection ActorMutableStateInspection
  private var lastSuccessOption: Option[OffsetDateTime] = None

  /**
    * Process the data asynchronously
    *
    * @return the number of elements processed
    */
  override protected def process(data: NonEmptyVector[(WorkflowId, OffsetDateTime)]): Future[Int] = instrumentedProcess {
    val heartbeatDateTime = OffsetDateTime.now()
    val processFuture = workflowStoreAccess.writeWorkflowHeartbeats(data, heartbeatDateTime)
    processFuture transform {
      // Track the `Try`, and then return the original `Try`. Similar to `andThen` but doesn't swallow exceptions.
      _ <| trackRepeatedFailures(heartbeatDateTime, data.length)
    }
  }

  override def receive = enhancedReceive.orElse(super.receive)
  override protected def weightFunction(command: (WorkflowId, OffsetDateTime)) = 1
  override protected def instrumentationPath = NonEmptyList.of("store", "heartbeat-writes")
  override protected def instrumentationPrefix = InstrumentationPrefixes.WorkflowPrefix
  override def commandToData(snd: ActorRef): PartialFunction[Any, (WorkflowId, OffsetDateTime)] = {
    case command: WorkflowStoreWriteHeartbeatCommand => (command.workflowId, command.submissionTime)
  }

  /*
  WARNING: Even though this is in an actor, the logic deals with instances of Future that could complete in _any_ order,
  and even call this method at the same time from different threads.

  We are expecting the underlying FSM to ensure that the call to this method does NOT occur in parallel, waiting for
  the call to `process` to complete.
   */
  private def trackRepeatedFailures(heartbeatDateTime: OffsetDateTime, workflowCount: Int)(processTry: Try[Int]): Unit = {
    processTry match {
      case Success(_) =>
        lastSuccessOption = Option(heartbeatDateTime)
      case Failure(_: Exception) =>
        /*
        If this is the first time processing heartbeats, then initialize the "last success" to when this heartbeat check
        began. This allows configuring the failure shutdown as low as zero seconds.
         */
        val lastSuccess = lastSuccessOption getOrElse {
          lastSuccessOption = Option(heartbeatDateTime)
          heartbeatDateTime
        }
        val now = OffsetDateTime.now()
        val failureJDuration = JDuration.between(lastSuccess, now)
        if (failureJDuration.toNanos >= failureShutdownDuration.toNanos) {
          val failureUnits = failureShutdownDuration.unit
          val failureLength = FiniteDuration(failureJDuration.toNanos, TimeUnit.NANOSECONDS).toUnit(failureUnits)
          log.error(String.format(
            "Shutting down %s as at least %s of heartbeat write errors have occurred between %s and %s (%s %s)",
            workflowHeartbeatConfig.cromwellId,
            failureShutdownDuration,
            lastSuccess,
            now,
            failureLength.toString,
            failureUnits.toString.toLowerCase
          ))
          terminator.beginCromwellShutdown(WorkflowStoreHeartbeatWriteActor.Shutdown)
        }
        ()
      case Failure(throwable) =>
        log.error(throwable, s"Shutting down due to ${throwable.getMessage}")
        terminator.beginCromwellShutdown(WorkflowStoreHeartbeatWriteActor.Shutdown)
        ()
    }
  }

}

object WorkflowStoreHeartbeatWriteActor {
  object Shutdown extends CoordinatedShutdown.Reason

  def props(
             workflowStoreAccess: WorkflowStoreAccess,
             workflowHeartbeatConfig: WorkflowHeartbeatConfig,
             terminator: CromwellTerminator,
             serviceRegistryActor: ActorRef
           ): Props =
    Props(
      WorkflowStoreHeartbeatWriteActor(
        workflowStoreAccess = workflowStoreAccess,
        workflowHeartbeatConfig = workflowHeartbeatConfig,
        terminator = terminator,
        serviceRegistryActor = serviceRegistryActor
      )).withDispatcher(EngineDispatcher)
}
