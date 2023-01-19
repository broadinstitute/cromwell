package cromwell.backend.google.pipelines.batch

import akka.actor.ActorRef
//import com.google.auth.Credentials
//import cromwell.backend.io.WorkflowPaths
//import cromwell.backend.wfs.WorkflowPathBuilder.workflowPaths
import cromwell.core.WorkflowOptions
//import cromwell.core.{Dispatcher, WorkflowOptions}

import scala.concurrent.Future
//import com.google.auth.Credentials

//import scala.concurrent.Future
//import cromwell.backend.google.pipelines.common.{PipelinesApiConfiguration, PipelinesApiInitializationActorParams, PipelinesApiRuntimeAttributes}
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
//import cromwell.core.WorkflowOptions
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
  protected val gcpBatchConfiguration: GcpBatchConfiguration = batchParams.batchConfiguration

  protected val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions

  //private lazy val ioEc = context.system.dispatchers.lookup(Dispatcher.IoDispatcher)

  /*
  private lazy val gcsCredentials: Future[Credentials] = gcpBatchConfiguration
    .papiAttributes
    .auths
    .gcs
    .retryCredentials(workflowOptions, List(StorageScopes
      .DEVSTORAGE_FULL_CONTROL))
*/

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    GcpBatchRuntimeAttributes
      .runtimeAttributesBuilder(gcpBatchConfiguration)

  override lazy val workflowPaths: Future[GcpBatchWorkflowPaths] = for {
    validatedPathBuilders <- pathBuilders
  } yield new GcpBatchWorkflowPaths(
      workflowDescriptor, gcpBatchConfiguration, validatedPathBuilders)


/*
  override lazy val initializationData: GcpBackendInitializationData = GcpBackendInitializationData(
    workflowPaths = GcpBatchWorkflowPaths,
    runtimeAttributesBuilder = runtimeAttributesBuilder,
    gcpBatchConfiguration = gcpBatchConfiguration
  )
  */

  override lazy val initializationData: Future[GcpBackendInitializationData] = for {
    batchWorkflowPaths <- workflowPaths
  } yield GcpBackendInitializationData(
    workflowPaths = batchWorkflowPaths,
    runtimeAttributesBuilder = runtimeAttributesBuilder,
    gcpBatchConfiguration = gcpBatchConfiguration
  )
  //add in gcs credentials if necessary

}

object GcpBatchInitializationActor {

}
