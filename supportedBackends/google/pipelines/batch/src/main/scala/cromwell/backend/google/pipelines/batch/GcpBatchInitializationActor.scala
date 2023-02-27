package cromwell.backend.google.pipelines.batch

import akka.actor.ActorRef
import com.google.api.services.lifesciences.v2beta.CloudLifeSciencesScopes
//import cromwell.backend.io.WorkflowPaths
//import cromwell.backend.wfs.WorkflowPathBuilder.workflowPaths
import cromwell.core.WorkflowOptions
//import cromwell.core.{Dispatcher, WorkflowOptions}

import scala.concurrent.Future
import com.google.auth.Credentials

//import scala.concurrent.Future
//import cromwell.backend.google.pipelines.common.{PipelinesApiConfiguration, PipelinesApiInitializationActorParams, PipelinesApiRuntimeAttributes}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.backend.standard.{StandardInitializationActor, StandardInitializationActorParams, StandardValidatedRuntimeAttributesBuilder}
//import cromwell.core.WorkflowOptions
import cromwell.core.io.AsyncIoActorClient
//import cromwell.core.io.{AsyncIoActorClient, IoCommandBuilder}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wom.graph.CommandCallNode
import com.google.api.services.storage.StorageScopes
import cromwell.filesystems.gcs.GoogleUtil._
import com.google.api.services.genomics.v2alpha1.GenomicsScopes

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
  protected val gcpBatchConfiguration: GcpBatchConfiguration = batchParams.batchConfiguration
  protected val workflowOptions: WorkflowOptions = workflowDescriptor.workflowOptions


  //private lazy val ioEc = context.system.dispatchers.lookup(Dispatcher.IoDispatcher)

  // Credentials object for the GCS API
  private lazy val gcsCredentials: Future[Credentials] = gcpBatchConfiguration
    .batchAttributes
    .auths
    .gcs
    .retryCredentials(workflowOptions, List(StorageScopes
      .DEVSTORAGE_FULL_CONTROL))

  // Credentials object for the Genomics API
  private lazy val genomicsCredentials: Future[Credentials] = gcpBatchConfiguration.batchAttributes.auths.genomics
                                                                                    .retryCredentials(workflowOptions, List(
                                                                                      CloudLifeSciencesScopes
                                                                                        .CLOUD_PLATFORM,
                                                                                      GenomicsScopes.GENOMICS
                                                                                    ))

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    GcpBatchRuntimeAttributes
      .runtimeAttributesBuilder(gcpBatchConfiguration)

  override lazy val workflowPaths: Future[GcpBatchWorkflowPaths] = for {
    gcsCred <- gcsCredentials
    genomicsCred <- genomicsCredentials
    validatedPathBuilders <- pathBuilders
  } yield new GcpBatchWorkflowPaths(
      workflowDescriptor, gcsCred, genomicsCred, gcpBatchConfiguration, validatedPathBuilders)


  override lazy val initializationData: Future[GcpBackendInitializationData] = for {
    batchWorkflowPaths <- workflowPaths
  } yield GcpBackendInitializationData(
    workflowPaths = batchWorkflowPaths,
    runtimeAttributesBuilder = runtimeAttributesBuilder,
    gcpBatchConfiguration = gcpBatchConfiguration
  )
  //add in gcs credentials if necessary

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    for {
      paths <- workflowPaths
      _ = publishWorkflowRoot(paths.workflowRoot.pathAsString)
      data <- initializationData
    } yield Option(data)
  }

  override lazy val ioCommandBuilder: GcsBatchCommandBuilder.type = GcsBatchCommandBuilder

}

object GcpBatchInitializationActor {

}
