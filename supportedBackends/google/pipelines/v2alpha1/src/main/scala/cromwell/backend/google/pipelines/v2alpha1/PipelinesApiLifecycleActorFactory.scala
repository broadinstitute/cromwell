package cromwell.backend.google.pipelines.v2alpha1

import akka.actor.ActorRef
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.common.{PipelinesApiBackendLifecycleActorFactory, PipelinesApiBackendSingletonActor, PipelinesApiConfiguration}
import cromwell.backend.google.pipelines.v2alpha1.api.request.RequestHandler
import cromwell.backend.standard.StandardAsyncExecutionActor

class PipelinesApiLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
  extends PipelinesApiBackendLifecycleActorFactory(name, configurationDescriptor) {
  val genomicsFactory = GenomicsFactory(googleConfig.applicationName, papiAttributes.auths.genomics, papiAttributes.endpointUrl)(papiAttributes.localizationConfiguration)
  override protected val jesConfiguration = new PipelinesApiConfiguration(configurationDescriptor, genomicsFactory, googleConfig, papiAttributes)
  override def requiredBackendSingletonActor(serviceRegistryActor: ActorRef) = {
    implicit val batchHandler = new RequestHandler(googleConfig.applicationName, papiAttributes.endpointUrl, papiAttributes)
    PipelinesApiBackendSingletonActor.props(
      jesConfiguration.papiAttributes.qps,
      jesConfiguration.papiAttributes.requestWorkers,
      serviceRegistryActor
    )
  }
  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[PipelinesApiAsyncBackendJobExecutionActor]
}
