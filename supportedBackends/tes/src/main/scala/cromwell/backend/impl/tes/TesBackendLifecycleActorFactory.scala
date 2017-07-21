package cromwell.backend.impl.tes

import akka.actor.ActorRef
import cromwell.backend._
import cromwell.backend.standard._
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import net.ceedubs.ficus.Ficus._
import wdl4s.wdl.WdlTaskCall

case class TesBackendLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {

  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[TesInitializationActor]

  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[TesAsyncBackendJobExecutionActor]

  override def jobIdKey: String = TesAsyncBackendJobExecutionActor.JobIdKey

  val tesConfiguration = new TesConfiguration(configurationDescriptor)

  override val jobExecutionTokenType: JobExecutionTokenType = {
    val concurrentJobLimit = configurationDescriptor.backendConfig.as[Option[Int]]("concurrent-job-limit")
    JobExecutionTokenType(name, concurrentJobLimit)
  }

  override def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[WdlTaskCall],
                                                 serviceRegistryActor: ActorRef): StandardInitializationActorParams = {
    TesInitializationActorParams(workflowDescriptor, calls, tesConfiguration, serviceRegistryActor)
  }
}
