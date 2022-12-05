package cromwell.backend.google.pipelines.batch

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.standard._


class GcpBatchBackendLifecycleActorFactory(name: String, override val configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {

  override def name: String = "batch"

  override def jobIdKey: String = "gcp_batch"

  override def asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[GcpBatchAsyncBackendJobExecutionActor]

  override def backendSingletonActorProps(serviceRegistryActor: ActorRef): Option[Props] = super
    .backendSingletonActorProps(serviceRegistryActor)

}

object GcpBatchBackendLifecycleActorFactory {

}
