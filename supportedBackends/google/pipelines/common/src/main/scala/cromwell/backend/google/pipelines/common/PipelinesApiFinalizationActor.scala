package cromwell.backend.google.pipelines.common

import akka.actor.ActorRef
import cromwell.backend._
import cromwell.backend.standard.{StandardFinalizationActor, StandardFinalizationActorParams}
import cromwell.core.CallOutputs
import cromwell.core.io.AsyncIoActorClient
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wom.graph.CommandCallNode

case class PipelinesApiFinalizationActorParams(
  workflowDescriptor: BackendWorkflowDescriptor,
  ioActor: ActorRef,
  calls: Set[CommandCallNode],
  jesConfiguration: PipelinesApiConfiguration,
  jobExecutionMap: JobExecutionMap,
  workflowOutputs: CallOutputs,
  initializationDataOption: Option[BackendInitializationData]
) extends StandardFinalizationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
}

class PipelinesApiFinalizationActor(val pipelinesParams: PipelinesApiFinalizationActorParams)
    extends StandardFinalizationActor(pipelinesParams)
    with AsyncIoActorClient {

  lazy val jesConfiguration: PipelinesApiConfiguration = pipelinesParams.jesConfiguration

  override lazy val ioCommandBuilder = GcsBatchCommandBuilder

  override def ioActor: ActorRef = pipelinesParams.ioActor
}
