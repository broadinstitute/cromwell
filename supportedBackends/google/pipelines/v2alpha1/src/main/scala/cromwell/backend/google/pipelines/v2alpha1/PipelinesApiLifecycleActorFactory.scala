package cromwell.backend.google.pipelines.v2alpha1

import akka.actor.ActorRef
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.google.pipelines.common.{PipelinesApiBackendLifecycleActorFactory, PipelinesApiBackendSingletonActor, PipelinesApiConfiguration}
import cromwell.backend.google.pipelines.v2alpha1.api.request.RequestHandler
import cromwell.backend.standard.StandardAsyncExecutionActor

class PipelinesApiLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
  extends PipelinesApiBackendLifecycleActorFactory(name, configurationDescriptor) {
  val genomicsFactory = GenomicsFactory(googleConfig.applicationName, jesAttributes.auths.genomics, jesAttributes.endpointUrl)
  override protected val jesConfiguration = PipelinesApiConfiguration(configurationDescriptor, genomicsFactory, googleConfig, jesAttributes)
  override def requiredBackendSingletonActor(serviceRegistryActor: ActorRef) = {
    implicit val batchHandler = new RequestHandler(googleConfig.applicationName, jesAttributes.endpointUrl)
    PipelinesApiBackendSingletonActor.props(jesConfiguration.qps, jesConfiguration.papiRequestWorkers, serviceRegistryActor)
  }
  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[PipelinesApiAsyncBackendJobExecutionActor]
}
