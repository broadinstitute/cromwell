package cromwell.backend.google.pipelines.v2alpha1

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.common.{
  PipelinesApiBackendLifecycleActorFactory,
  PipelinesApiBackendSingletonActor,
  PipelinesApiConfiguration
}
import cromwell.backend.google.pipelines.v2alpha1.api.request.RequestHandler
import cromwell.backend.standard.StandardAsyncExecutionActor

class PipelinesApiLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
    extends PipelinesApiBackendLifecycleActorFactory(name, configurationDescriptor) {
  private val genomicsFactory =
    GenomicsFactory(
      googleConfig.applicationName,
      papiAttributes.auths.genomics,
      papiAttributes.endpointUrl
    )(papiAttributes.gcsTransferConfiguration)
  override protected val jesConfiguration =
    new PipelinesApiConfiguration(configurationDescriptor, genomicsFactory, googleConfig, papiAttributes)
  override def requiredBackendSingletonActor(serviceRegistryActor: ActorRef): Props = {
    implicit val batchHandler: RequestHandler =
      new RequestHandler(
        googleConfig.applicationName,
        papiAttributes.endpointUrl,
        papiAttributes.batchRequestTimeoutConfiguration
      )
    PipelinesApiBackendSingletonActor.props(
      jesConfiguration.papiAttributes.qps,
      jesConfiguration.papiAttributes.requestWorkers,
      serviceRegistryActor
    )
  }
  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[PipelinesApiAsyncBackendJobExecutionActor]
}
