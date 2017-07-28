package cromwell.backend.impl.jes

import akka.actor.ActorRef
import cromwell.backend._
import cromwell.backend.impl.jes.callcaching.{JesBackendCacheHitCopyingActor, JesBackendFileHashingActor}
import cromwell.backend.standard._
import cromwell.backend.standard.callcaching.{StandardCacheHitCopyingActor, StandardFileHashingActor}
import cromwell.core.CallOutputs
import wdl4s.wdl.WdlTaskCall

import scala.util.{Success, Try}
import cromwell.backend.impl.jes.JesBackendLifecycleActorFactory._

case class JesBackendLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {

  override lazy val initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[JesInitializationActor]

  override lazy val asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[JesAsyncBackendJobExecutionActor]

  override lazy val finalizationActorClassOption: Option[Class[_ <: StandardFinalizationActor]] =
    Option(classOf[JesFinalizationActor])

  override lazy val jobIdKey: String = JesAsyncBackendJobExecutionActor.JesOperationIdKey

  val jesConfiguration = new JesConfiguration(configurationDescriptor)

  override def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[WdlTaskCall],
                                                 serviceRegistryActor: ActorRef, restart: Boolean): StandardInitializationActorParams = {
    JesInitializationActorParams(workflowDescriptor, ioActor, calls, jesConfiguration, serviceRegistryActor, restart)
  }

  override def workflowFinalizationActorParams(workflowDescriptor: BackendWorkflowDescriptor, ioActor: ActorRef, calls: Set[WdlTaskCall],
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
    Option(classOf[JesBackendCacheHitCopyingActor])
  }

  override def backendSingletonActorProps = Option(JesBackendSingletonActor.props(jesConfiguration.qps))

  override lazy val fileHashingActorClassOption: Option[Class[_ <: StandardFileHashingActor]] = Option(classOf[JesBackendFileHashingActor])
  
  override def dockerHashCredentials(initializationData: Option[BackendInitializationData]) = {
    Try(BackendInitializationData.as[JesBackendInitializationData](initializationData)) match {
      case Success(jesData) =>
        val maybeDockerHubCredentials = jesData.jesConfiguration.dockerCredentials
        val googleCredentials = Option(jesData.gcsCredentials)
        List(maybeDockerHubCredentials, googleCredentials).flatten
      case _ => List.empty[Any]
    }
  }
  override val requestedKeyValueStoreKeys: Seq[String] = Seq(preemptionCountKey, unexpectedRetryCountKey)
}

object JesBackendLifecycleActorFactory {
  val preemptionCountKey = "PreemptionCount"
  val unexpectedRetryCountKey = "UnexpectedRetryCount"
}
