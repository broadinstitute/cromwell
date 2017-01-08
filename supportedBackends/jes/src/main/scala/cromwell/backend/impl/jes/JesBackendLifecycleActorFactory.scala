package cromwell.backend.impl.jes

import akka.actor.ActorRef
import cromwell.backend._
import cromwell.backend.callcaching.FileHashingActor.FileHashingFunction
import cromwell.backend.impl.jes.callcaching.JesBackendFileHashing
import cromwell.backend.standard._
import cromwell.core.CallOutputs
import wdl4s.TaskCall
import wdl4s.expression.WdlStandardLibraryFunctions

case class JesBackendLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
  extends StandardLifecycleActorFactory {

  override def initializationActorClass: Class[_ <: StandardInitializationActor] = classOf[JesInitializationActor]

  override def asyncExecutionActorClass: Class[_ <: StandardAsyncExecutionActor] =
    classOf[JesAsyncBackendJobExecutionActor]

  override def standardCacheHitCopyingActorOption: Option[Class[_ <: StandardCacheHitCopyingActor]] =
    Option(classOf[JesCacheHitCopyingActor])

  override def finalizationActorClassOption: Option[Class[_ <: StandardFinalizationActor]] =
    Option(classOf[JesFinalizationActor])

  override def jobIdKey: String = JesAsyncBackendJobExecutionActor.JesOperationIdKey

  val jesConfiguration = new JesConfiguration(configurationDescriptor)

  override def workflowInitializationActorParams(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[TaskCall],
                                                 serviceRegistryActor: ActorRef): StandardInitializationActorParams = {
    JesInitializationActorParams(workflowDescriptor, calls, jesConfiguration, serviceRegistryActor)
  }

  override def workflowFinalizationActorParams(workflowDescriptor: BackendWorkflowDescriptor, calls: Set[TaskCall],
                                              jobExecutionMap: JobExecutionMap, workflowOutputs: CallOutputs,
                                              initializationDataOption: Option[BackendInitializationData]):
  StandardFinalizationActorParams = {
    // The `JesInitializationActor` will only return a non-`Empty` `JesBackendInitializationData` from a successful `beforeAll`
    // invocation.  HOWEVER, the finalization actor is created regardless of whether workflow initialization was successful
    // or not.  So the finalization actor must be able to handle an empty `JesBackendInitializationData` option, and there is no
    // `.get` on the initialization data as there is with the execution or cache hit copying actor methods.
    JesFinalizationActorParams(workflowDescriptor, calls, jesConfiguration, jobExecutionMap, workflowOutputs,
      initializationDataOption)
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationDataOption: Option[BackendInitializationData]
                                          ): WdlStandardLibraryFunctions = {
    val initializationData = BackendInitializationData.
      as[JesBackendInitializationData](initializationDataOption)

    val jesCallPaths = initializationData.workflowPaths.toJobPaths(jobKey, workflowDescriptor)
    new JesExpressionFunctions(List(jesCallPaths.gcsPathBuilder), jesCallPaths.callContext)
  }

  override def backendSingletonActorProps = Option(JesBackendSingletonActor.props(jesConfiguration.qps))

  override lazy val fileHashingFunction: Option[FileHashingFunction] = Option(FileHashingFunction(JesBackendFileHashing.getCrc32c))
}
