package cromwell.backend.impl.tes

import akka.actor.ActorRef
import cromwell.backend._
import cromwell.backend.standard._
import cromwell.backend.standard.callcaching.StandardFileHashingActor
import wom.graph.CommandCallNode

case class TesBackendLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
    extends StandardLifecycleActorFactory
    with AzurePlatform {

  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[TesInitializationActor]

  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[TesAsyncBackendJobExecutionActor]

  override def jobIdKey: String = TesAsyncBackendJobExecutionActor.JobIdKey

  val tesConfiguration = new TesConfiguration(configurationDescriptor)

  override def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor,
                                                 ioActor: ActorRef,
                                                 calls: Set[CommandCallNode],
                                                 serviceRegistryActor: ActorRef,
                                                 restarting: Boolean
  ): StandardInitializationActorParams =
    TesInitializationActorParams(workflowDescriptor, calls, tesConfiguration, serviceRegistryActor)

  // Use our custom hashing actor so that local-style paths (/mnt/efs/…, /some/path)
  // are hashed with sibling-MD5 support instead of the default async IO hash.
  override lazy val fileHashingActorClassOption: Option[Class[_ <: StandardFileHashingActor]] =
    Option(classOf[TesBackendFileHashingActor])
}
