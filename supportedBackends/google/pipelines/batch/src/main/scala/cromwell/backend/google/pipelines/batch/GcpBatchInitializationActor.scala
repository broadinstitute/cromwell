package cromwell.backend.google.pipelines.batch

import akka.actor.ActorRef
//import cromwell.backend.google.pipelines.common.{PipelinesApiConfiguration, PipelinesApiInitializationActorParams, PipelinesApiRuntimeAttributes}
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
import cromwell.core.WorkflowOptions
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
class GcpBatchInitializationActor(batchParams: GcpBatchInitializationActorParams) extends StandardInitializationActor(batchParams) with AsyncIoActorClient {
  override lazy val ioActor: ActorRef = batchParams.ioActor

  override lazy val ioCommandBuilder: GcsBatchCommandBuilder.type = GcsBatchCommandBuilder
  protected val batchConfiguration: GcpBatchConfiguration = batchParams.batchConfiguration

  protected val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    GcpBatchRuntimeAttributes
      .runtimeAttributesBuilder(batchConfiguration)
}

object GcpBatchInitializationActor {

}
