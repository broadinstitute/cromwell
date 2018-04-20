package cromwell.backend.google.pipelines.common

import akka.actor.ActorRef
import cromwell.backend._
import cromwell.backend.google.pipelines.common.PipelinesApiBackendLifecycleActorFactory._
import cromwell.backend.google.pipelines.common.callcaching.{PipelinesApiBackendCacheHitCopyingActor, PipelinesApiBackendFileHashingActor}
import cromwell.backend.standard._
import cromwell.backend.standard.callcaching.{StandardCacheHitCopyingActor, StandardFileHashingActor}
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.core.CallOutputs
import wom.graph.CommandCallNode

import scala.util.{Success, Try}

abstract class PipelinesApiBackendLifecycleActorFactory(override val name: String, override val configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {

  override val requestedKeyValueStoreKeys: Seq[String] = Seq(preemptionCountKey, unexpectedRetryCountKey)

  protected val googleConfig = GoogleConfiguration(configurationDescriptor.globalConfig)
  protected val jesAttributes = PipelinesApiAttributes(googleConfig, configurationDescriptor.backendConfig)

  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[JesInitializationActor]

  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[PipelinesApiAsyncBackendJobExecutionActor]

  override lazy val finalizationActorClassOption: Option[Class[_ <: StandardFinalizationActor]] =
    Option(classOf[JesFinalizationActor])

  override lazy val jobIdKey: String = PipelinesApiAsyncBackendJobExecutionActor.JesOperationIdKey

  protected val jesConfiguration: PipelinesApiConfiguration

  override def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[CommandCallNode],
                                                 serviceRegistryActor: ActorRef, restart: Boolean): StandardInitializationActorParams = {
    JesInitializationActorParams(workflowDescriptor, ioActor, calls, jesConfiguration, serviceRegistryActor, restart)
  }

  override def workflowFinalizationActorParams(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[CommandCallNode],
                                               jobExecutionMap: JobExecutionMap, workflowOutputs: CallOutputs,
                                               initializationDataOption: Option[BackendInitializationData]):
  StandardFinalizationActorParams = {
    // The `JesInitializationActor` will only return a non-`Empty` `JesBackendInitializationData` from a successful `beforeAll`
    // invocation.  HOWEVER, the finalization actor is created regardless of whether workflow initialization was successful
    // or not.  So the finalization actor must be able to handle an empty `JesBackendInitializationData` option, and there is no
    // `.get` on the initialization data as there is with the execution or cache hit copying actor methods.
    JesFinalizationActorParams(workflowDescriptor, ioActor, calls, jesConfiguration, jobExecutionMap, workflowOutputs,
      initializationDataOption)
  }

  override lazy val cacheHitCopyingActorClassOption: Option[Class[_ <: StandardCacheHitCopyingActor]] = {
    Option(classOf[PipelinesApiBackendCacheHitCopyingActor])
  }

  override lazy val fileHashingActorClassOption: Option[Class[_ <: StandardFileHashingActor]] = Option(classOf[PipelinesApiBackendFileHashingActor])

  override def dockerHashCredentials(initializationData: Option[BackendInitializationData]) = {
    Try(BackendInitializationData.as[PipelinesApiBackendInitializationData](initializationData)) match {
      case Success(jesData) =>
        val maybeDockerHubCredentials = jesData.jesConfiguration.dockerCredentials
        val googleCredentials = Option(jesData.gcsCredentials)
        List(maybeDockerHubCredentials, googleCredentials).flatten
      case _ => List.empty[Any]
    }
  }
}

object PipelinesApiBackendLifecycleActorFactory {
  val preemptionCountKey = "PreemptionCount"
  val unexpectedRetryCountKey = "UnexpectedRetryCount"
}
