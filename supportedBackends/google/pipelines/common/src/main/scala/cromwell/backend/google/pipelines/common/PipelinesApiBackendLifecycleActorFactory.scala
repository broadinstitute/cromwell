package cromwell.backend.google.pipelines.common

import akka.actor.{ActorRef, Props}
import cromwell.backend._
import cromwell.backend.google.pipelines.common.PipelinesApiBackendLifecycleActorFactory._
import cromwell.backend.google.pipelines.common.authentication.PipelinesApiDockerCredentials
import cromwell.backend.google.pipelines.common.callcaching.{PipelinesApiBackendCacheHitCopyingActor, PipelinesApiBackendFileHashingActor}
import cromwell.backend.standard._
import cromwell.backend.standard.callcaching.{StandardCacheHitCopyingActor, StandardFileHashingActor}
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.{CallOutputs, DockerCredentials}
import wom.graph.CommandCallNode

import scala.util.{Success, Try}

abstract class PipelinesApiBackendLifecycleActorFactory(override val name: String, override val configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {

  // Abstract members
  protected def requiredBackendSingletonActor(serviceRegistryActor: ActorRef): Props
  protected val jesConfiguration: PipelinesApiConfiguration

  override val requestedKeyValueStoreKeys: Seq[String] = Seq(preemptionCountKey, unexpectedRetryCountKey)

  protected val googleConfig = GoogleConfiguration(configurationDescriptor.globalConfig)

  protected val papiAttributes = PipelinesApiConfigurationAttributes(googleConfig, configurationDescriptor.backendConfig)

  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[PipelinesApiInitializationActor]

  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[PipelinesApiAsyncBackendJobExecutionActor]
  override lazy val finalizationActorClassOption: Option[Class[_ <: StandardFinalizationActor]] =
    Option(classOf[PipelinesApiFinalizationActor])
  override lazy val jobIdKey: String = PipelinesApiAsyncBackendJobExecutionActor.JesOperationIdKey

  override def backendSingletonActorProps(serviceRegistryActor: ActorRef) = Option(requiredBackendSingletonActor(serviceRegistryActor))

  override def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[CommandCallNode],
                                                 serviceRegistryActor: ActorRef, restart: Boolean): StandardInitializationActorParams = {
    PipelinesApiInitializationActorParams(workflowDescriptor, ioActor, calls, jesConfiguration, serviceRegistryActor, restart)
  }

  override def workflowFinalizationActorParams(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[CommandCallNode],
                                               jobExecutionMap: JobExecutionMap, workflowOutputs: CallOutputs,
                                               initializationDataOption: Option[BackendInitializationData]):
  StandardFinalizationActorParams = {
    // The `PipelinesApiInitializationActor` will only return a non-`Empty` `PipelinesApiBackendInitializationData` from a successful `beforeAll`
    // invocation.  HOWEVER, the finalization actor is created regardless of whether workflow initialization was successful
    // or not.  So the finalization actor must be able to handle an empty `PipelinesApiBackendInitializationData` option, and there is no
    // `.get` on the initialization data as there is with the execution or cache hit copying actor methods.
    PipelinesApiFinalizationActorParams(workflowDescriptor, ioActor, calls, jesConfiguration, jobExecutionMap, workflowOutputs,
      initializationDataOption)
  }

  override lazy val cacheHitCopyingActorClassOption: Option[Class[_ <: StandardCacheHitCopyingActor]] = {
    Option(classOf[PipelinesApiBackendCacheHitCopyingActor])
  }

  override lazy val fileHashingActorClassOption: Option[Class[_ <: StandardFileHashingActor]] = Option(classOf[PipelinesApiBackendFileHashingActor])

  override def dockerHashCredentials(workflowDescriptor: BackendWorkflowDescriptor, initializationData: Option[BackendInitializationData]) = {
    Try(BackendInitializationData.as[PipelinesApiBackendInitializationData](initializationData)) match {
      case Success(papiData) =>
        val tokenFromWorkflowOptions = workflowDescriptor.workflowOptions.get(GoogleAuthMode.DockerCredentialsTokenKey).toOption
        val effectiveToken = tokenFromWorkflowOptions.orElse(papiData.papiConfiguration.dockerCredentials map { _.token })

        val dockerCredentials: Option[PipelinesApiDockerCredentials] = effectiveToken map { token =>
          // These credentials are being returned for hashing and all that matters in this context is the token
          // so just `None` the auth and key.
          val baseDockerCredentials = new DockerCredentials(token = token, authName = None, keyName = None)
          PipelinesApiDockerCredentials.apply(baseDockerCredentials, googleConfig)
        }
        val googleCredentials = Option(papiData.gcsCredentials)
        List(dockerCredentials, googleCredentials).flatten
      case _ => List.empty[Any]
    }
  }
}

object PipelinesApiBackendLifecycleActorFactory {
  val preemptionCountKey = "PreemptionCount"
  val unexpectedRetryCountKey = "UnexpectedRetryCount"
}
