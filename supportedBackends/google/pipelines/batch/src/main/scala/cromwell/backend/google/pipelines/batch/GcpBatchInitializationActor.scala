package cromwell.backend.google.pipelines.batch

import akka.actor.ActorRef
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams}
import cromwell.core.io.AsyncIoActorClient
//import cromwell.core.io.{AsyncIoActorClient, IoCommandBuilder}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wom.graph.CommandCallNode

case class GcpBatchInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  ioActor: ActorRef,
  calls: Set[CommandCallNode],
  batchConfiguration: GcpBatchConfiguration,
  serviceRegistryActor: ActorRef,
  restarting: Boolean
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = batchConfiguration.configurationDescriptor

}
class GcpBatchInitializationActor(gcpBatchParams: GcpBatchInitializationActorParams) extends StandardInitializationActor(gcpBatchParams) with AsyncIoActorClient {
  override lazy val ioActor: ActorRef = gcpBatchParams.ioActor

  override lazy val ioCommandBuilder: GcsBatchCommandBuilder.type = GcsBatchCommandBuilder
}

object GcpBatchInitializationActor {

}
