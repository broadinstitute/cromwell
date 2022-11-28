package cromwell.backend.google.pipelines.batch

//import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.standard._

class GcpBatchBackendLifecycleActorFactory(name: String, override val configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {

  // Abstract members
 // protected def requiredBackendSingletonActor(serviceRegistryActor: ActorRef): Props
  //override def backendSingletonActorProps(serviceRegistryActor: ActorRef): Option[Props] = Option(requiredBackendSingletonActor(serviceRegistryActor))

  //override def configurationDescriptor: BackendConfigurationDescriptor = GcpBatchConfiguration(configurationDescriptor)
  //val gcpBatchConfig = new GcpBatchConfiguration(configurationDescriptor)

  override def jobIdKey: String = "gcp_batch"

  //override def name: String = "batch"

   //override def configurationDescriptor: BackendConfigurationDescriptor = new GcpBatchConfiguration(configurationDescriptor)

  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[GcpBatchAsyncBackendJobExecutionActor]


  override def name: String = "batch"
}

object GcpBatchBackendLifecycleActorFactory {

}
