package cromwell.backend.google.pipelines.batch

import akka.actor.{ActorRef, Props}
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend.BackendWorkflowDescriptor
import cromwell.backend.{BackendInitializationData, JobExecutionMap}
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.CallOutputs
import wom.graph.CommandCallNode
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.standard._


class GcpBatchBackendLifecycleActorFactory(name: String, override val configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {

  override def name: String = "batch"

  override def jobIdKey: String = "__gcp_batch"
  protected val googleConfig: GoogleConfiguration = GoogleConfiguration(configurationDescriptor.globalConfig)

  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[GcpBatchInitializationActor]

  override def asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[GcpBatchAsyncBackendJobExecutionActor]

  override lazy val finalizationActorClassOption: Option[Class[_ <: StandardFinalizationActor]] =
    Option(classOf[GcpBatchFinalizationActor])


  val batchConfiguration = new GcpBatchConfiguration(configurationDescriptor, googleConfig)


  override def workflowInitializationActorParams(

                                                  workflowDescriptor: BackendWorkflowDescriptor,
                                                  ioActor: ActorRef,
                                                  calls: Set[CommandCallNode],
                                                  serviceRegistryActor: ActorRef,
                                                  restart: Boolean): StandardInitializationActorParams = {
    GcpBatchInitializationActorParams(workflowDescriptor, ioActor , calls, batchConfiguration, serviceRegistryActor, restart)
  }

  override def workflowFinalizationActorParams(
                                           workflowDescriptor: BackendWorkflowDescriptor,
                                           ioActor: ActorRef,
                                           //batchConfiguration: GcpBatchConfiguration,
                                           calls: Set[CommandCallNode],
                                           jobExecutionMap: JobExecutionMap,
                                           workflowOutputs: CallOutputs,
                                           initializationDataOption: Option[BackendInitializationData]): StandardFinalizationActorParams = {
    GcpBatchFinalizationActorParams(workflowDescriptor, ioActor, batchConfiguration, calls, jobExecutionMap, workflowOutputs, initializationDataOption)
  }

  override def backendSingletonActorProps(serviceRegistryActor: ActorRef): Option[Props] = {
    Option(GcpBatchBackendSingletonActor.props("gcp-batch"))
  }

  //override def backendSingletonActorProps(serviceRegistryActor: ActorRef): Option[Props] = super
  //  .backendSingletonActorProps(serviceRegistryActor)

}

object GcpBatchBackendLifecycleActorFactory extends StrictLogging {


}
