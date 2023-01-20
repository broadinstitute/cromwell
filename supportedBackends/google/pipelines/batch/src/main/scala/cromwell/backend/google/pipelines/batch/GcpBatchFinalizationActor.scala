package cromwell.backend.google.pipelines.batch

import akka.actor.ActorRef
import cromwell.backend._
//import cromwell.backend.google.pipelines.common.{PipelinesApiConfiguration, PipelinesApiFinalizationActorParams}
import cromwell.backend.standard.{StandardFinalizationActor, StandardFinalizationActorParams}
import cromwell.core.CallOutputs
import cromwell.core.io.AsyncIoActorClient
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wom.graph.CommandCallNode
case class GcpBatchFinalizationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  ioActor: ActorRef,
  batchConfiguration: GcpBatchConfiguration,
  calls: Set[CommandCallNode],
  jobExecutionMap: JobExecutionMap,
  workflowOutputs: CallOutputs,
  initializationDataOption: Option[BackendInitializationData]
) extends StandardFinalizationActorParams {
  //override def initializationDataOption: Option[BackendInitializationData] = ???
  override def configurationDescriptor: BackendConfigurationDescriptor = batchConfiguration.configurationDescriptor
}

class GcpBatchFinalizationActor(val batchParams: GcpBatchFinalizationActorParams) extends StandardFinalizationActor(batchParams) with AsyncIoActorClient {

  //lazy val batchConfiguration: GcpBatchConfiguration = gcpBatchParams.batchConfiguration
  lazy val batchConfiguration: GcpBatchConfiguration = batchParams.batchConfiguration

  override lazy val ioCommandBuilder = GcsBatchCommandBuilder
  override def ioActor: ActorRef = batchParams.ioActor
}

object GcpBatchFinalizationActor {

}