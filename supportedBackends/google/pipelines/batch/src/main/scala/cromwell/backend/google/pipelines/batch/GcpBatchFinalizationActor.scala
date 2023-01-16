package cromwell.backend.google.pipelines.batch

import akka.actor.ActorRef
import cromwell.backend._
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

class GcpBatchFinalizationActor(val gcpBatchParams: GcpBatchFinalizationActorParams) extends StandardFinalizationActor(gcpBatchParams) with AsyncIoActorClient {

  lazy val batchConfiguration: GcpBatchConfiguration = gcpBatchParams.batchConfiguration

  override lazy val ioCommandBuilder = GcsBatchCommandBuilder
  override def ioActor: ActorRef = gcpBatchParams.ioActor
}

object GcpBatchFinalizationActor {

}