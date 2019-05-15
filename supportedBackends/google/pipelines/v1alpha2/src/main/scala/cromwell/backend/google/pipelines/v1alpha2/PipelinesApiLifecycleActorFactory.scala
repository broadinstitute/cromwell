package cromwell.backend.google.pipelines.v1alpha2

import akka.actor.ActorRef
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.common.{PipelinesApiBackendLifecycleActorFactory, PipelinesApiBackendSingletonActor, PipelinesApiConfiguration}
import cromwell.backend.google.pipelines.v1alpha2.api.request.RequestHandler
import cromwell.backend.standard.{StandardFinalizationActor, StandardInitializationActor}

class PipelinesApiLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
  extends PipelinesApiBackendLifecycleActorFactory(name, configurationDescriptor) {
  private val genomicsFactory = GenomicsFactory(googleConfig.applicationName, papiAttributes.auths.genomics, papiAttributes.endpointUrl)
  override protected val jesConfiguration = new PipelinesApiConfiguration(configurationDescriptor, genomicsFactory, googleConfig, papiAttributes)
  override def requiredBackendSingletonActor(serviceRegistryActor: ActorRef) = {
    implicit val batchHandler = new RequestHandler(googleConfig.applicationName, papiAttributes.endpointUrl)
    PipelinesApiBackendSingletonActor.props(
      jesConfiguration.papiAttributes.qps,
      jesConfiguration.papiAttributes.requestWorkers,
      serviceRegistryActor
    )
  }

  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] =
    classOf[PipelinesApiInitializationActor]

  override lazy val finalizationActorClassOption: Option[Class[_ <: StandardFinalizationActor]] =
    Option(classOf[PipelinesApiFinalizationActor])
}
