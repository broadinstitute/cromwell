package cromwell.backend.impl.bcs

import akka.actor.ActorRef
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor}
import cromwell.backend.standard._
import cromwell.backend.BackendInitializationData
import cromwell.backend.impl.bcs.callcaching.BcsBackendCacheHitCopyingActor
import cromwell.backend.standard.callcaching.StandardCacheHitCopyingActor
import wom.graph.CommandCallNode

import scala.util.{Success, Try}


final case class BcsBackendLifecycleActorFactory(val name: String, val configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {
  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[BcsInitializationActor]
  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] = classOf[BcsAsyncBackendJobExecutionActor]

  override def jobIdKey: String = BcsAsyncBackendJobExecutionActor.JobIdKey

  val bcsConfiguration = new BcsConfiguration(configurationDescriptor)

  override def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[CommandCallNode], serviceRegistryActor: ActorRef, restarting: Boolean): StandardInitializationActorParams = {
    BcsInitializationActorParams(workflowDescriptor, calls, bcsConfiguration, serviceRegistryActor)
  }

  override lazy val cacheHitCopyingActorClassOption: Option[Class[_ <: StandardCacheHitCopyingActor]] = {
    Option(classOf[BcsBackendCacheHitCopyingActor])
  }

  override def dockerHashCredentials(workflowDescriptor: BackendWorkflowDescriptor, initializationData: Option[BackendInitializationData]) = {
    Try(BackendInitializationData.as[BcsBackendInitializationData](initializationData)) match {
      case Success(bcsData) =>
        bcsData.bcsConfiguration.dockerHashEndpoint match {
          case Some(endpoint) => List(bcsData.bcsConfiguration.dockerCredentials, Option(endpoint)).flatten
          case None => List(bcsData.bcsConfiguration.dockerCredentials).flatten
        }
      case _ => List.empty[Any]
    }
  }  
}
