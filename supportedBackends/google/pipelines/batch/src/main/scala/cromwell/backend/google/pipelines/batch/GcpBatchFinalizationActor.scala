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
  calls: Set[CommandCallNode],
  jobExecutionMap: JobExecutionMap,
  workflowOutputs: CallOutputs
) extends StandardFinalizationActorParams {
  override def initializationDataOption: Option[BackendInitializationData] = ???
  override def configurationDescriptor: BackendConfigurationDescriptor = ???
}

class GcpBatchFinalizationActor(val gcpBatchParams: GcpBatchFinalizationActorParams) extends StandardFinalizationActor(gcpBatchParams) with AsyncIoActorClient {

  //lazy val jesConfiguration: GcpBatchConfiguration = pipelinesParams.jesConfiguration

  override lazy val ioCommandBuilder = GcsBatchCommandBuilder
  override def ioActor: ActorRef = gcpBatchParams.ioActor
}