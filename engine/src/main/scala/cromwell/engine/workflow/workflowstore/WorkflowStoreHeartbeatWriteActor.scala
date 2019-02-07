package cromwell.engine.workflow.workflowstore

import java.time.OffsetDateTime

import akka.actor.{ActorRef, Props}
import cats.data.{NonEmptyList, NonEmptyVector}
import cromwell.core.Dispatcher.EngineDispatcher
import cromwell.core.WorkflowId
import cromwell.core.instrumentation.InstrumentationPrefixes
import cromwell.engine.workflow.workflowstore.WorkflowStoreActor.WorkflowStoreWriteHeartbeatCommand
import cromwell.services.EnhancedBatchActor

import scala.concurrent.Future

case class WorkflowStoreHeartbeatWriteActor(workflowStoreAccess: WorkflowStoreAccess,
                                            workflowHeartbeatConfig: WorkflowHeartbeatConfig,
                                            override val serviceRegistryActor: ActorRef)

  extends EnhancedBatchActor[(WorkflowId, OffsetDateTime)](
    flushRate = workflowHeartbeatConfig.heartbeatInterval,
    batchSize = workflowHeartbeatConfig.writeBatchSize) {

  override val threshold = workflowHeartbeatConfig.writeThreshold

  /**
    * Process the data asynchronously
    *
    * @return the number of elements processed
    */
  override protected def process(data: NonEmptyVector[(WorkflowId, OffsetDateTime)]): Future[Int] = instrumentedProcess {
    workflowStoreAccess.writeWorkflowHeartbeats(data)
  }

  override def receive = enhancedReceive.orElse(super.receive)
  override protected def weightFunction(command: (WorkflowId, OffsetDateTime)) = 1
  override protected def instrumentationPath = NonEmptyList.of("store", "heartbeat-writes")
  override protected def instrumentationPrefix = InstrumentationPrefixes.WorkflowPrefix
  override def commandToData(snd: ActorRef): PartialFunction[Any, (WorkflowId, OffsetDateTime)] = {
    case command: WorkflowStoreWriteHeartbeatCommand => (command.workflowId, command.submissionTime)
  }
}

object WorkflowStoreHeartbeatWriteActor {
  def props(
             workflowStoreAccess: WorkflowStoreAccess,
             workflowHeartbeatConfig: WorkflowHeartbeatConfig,
             serviceRegistryActor: ActorRef
           ): Props =
    Props(
      WorkflowStoreHeartbeatWriteActor(
        workflowStoreAccess = workflowStoreAccess,
        workflowHeartbeatConfig = workflowHeartbeatConfig,
        serviceRegistryActor = serviceRegistryActor
      )).withDispatcher(EngineDispatcher)
}
