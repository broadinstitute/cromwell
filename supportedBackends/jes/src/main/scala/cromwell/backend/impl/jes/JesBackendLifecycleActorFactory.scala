package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.{ActorRef, Props}
import com.typesafe.config.Config
import cromwell.backend._
import cromwell.backend.callcaching.FileHashingActor.FileHashingFunction
import cromwell.backend.impl.jes.callcaching.JesBackendFileHashing
import cromwell.backend.validation.RuntimeAttributesKeys
import cromwell.core.CallOutputs
import cromwell.core.Dispatcher.BackendDispatcher
import wdl4s.TaskCall
import wdl4s.expression.WdlStandardLibraryFunctions


case class JesBackendLifecycleActorFactory(name: String, configurationDescriptor: BackendConfigurationDescriptor)
  extends BackendLifecycleActorFactory {
  import JesBackendLifecycleActorFactory._

  val jesConfiguration = new JesConfiguration(configurationDescriptor)

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Set[TaskCall],
                                                serviceRegistryActor: ActorRef): Option[Props] = {
    Option(JesInitializationActor.props(workflowDescriptor, calls, jesConfiguration, serviceRegistryActor).withDispatcher(BackendDispatcher))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef,
                                      backendSingletonActor: Option[ActorRef]): Props = {
    // The `JesInitializationActor` will only return a non-`Empty` `JesBackendInitializationData` from a successful `beforeAll`
    // invocation, so the `get` here is safe.
    JesJobExecutionActor.props(jobDescriptor, jesConfiguration, initializationData.toJes.get, serviceRegistryActor, backendSingletonActor).withDispatcher(BackendDispatcher)
  }

  override def cacheHitCopyingActorProps = Option(cacheHitCopyingActorInner _)

  def cacheHitCopyingActorInner(jobDescriptor: BackendJobDescriptor,
                                         initializationData: Option[BackendInitializationData],
                                         serviceRegistryActor: ActorRef): Props = {
    // The `JesInitializationActor` will only return a non-`Empty` `JesBackendInitializationData` from a successful `beforeAll`
    // invocation, so the `get` here is safe.
    JesCacheHitCopyingActor.props(jobDescriptor, jesConfiguration, initializationData.toJes.get, serviceRegistryActor).withDispatcher(BackendDispatcher)
  }

  override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                              calls: Set[TaskCall],
                                              jobExecutionMap: JobExecutionMap,
                                              workflowOutputs: CallOutputs,
                                              initializationData: Option[BackendInitializationData]) = {
    // The `JesInitializationActor` will only return a non-`Empty` `JesBackendInitializationData` from a successful `beforeAll`
    // invocation.  HOWEVER, the finalization actor is created regardless of whether workflow initialization was successful
    // or not.  So the finalization actor must be able to handle an empty `JesBackendInitializationData` option, and there is no
    // `.get` on the initialization data as there is with the execution or cache hit copying actor methods.
    Option(JesFinalizationActor.props(workflowDescriptor, calls, jesConfiguration, jobExecutionMap, workflowOutputs, initializationData.toJes).withDispatcher(BackendDispatcher))
  }

  override def runtimeAttributeDefinitions(initializationDataOption: Option[BackendInitializationData]) = staticRuntimeAttributeDefinitions

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {

    val jesCallPaths = initializationData.toJes.get.workflowPaths.toJobPaths(jobKey)
    new JesExpressionFunctions(List(jesCallPaths.gcsPathBuilder), jesCallPaths.callContext)
  }

  override def getExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config,
                                    initializationData: Option[BackendInitializationData]): Path = {
    initializationData.toJes.get.workflowPaths.executionRoot
  }

  override def getWorkflowExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config,
                                    initializationData: Option[BackendInitializationData]): Path = {
    initializationData.toJes.get.workflowPaths.workflowRoot
  }

  override def backendSingletonActorProps = Option(JesBackendSingletonActor.props(jesConfiguration.qps))

  override lazy val fileHashingFunction: Option[FileHashingFunction] = Option(FileHashingFunction(JesBackendFileHashing.getCrc32c))
}

object JesBackendLifecycleActorFactory {
  implicit class Jessify(val genericInitializationData: Option[BackendInitializationData]) {
    // This leaves the result in an `Option` as finalization will be called even if initialization has failed, and if
    // initialization fails there won't be any initialization data.  The various `.get`s that occur below are in instances
    // where the workflow has successfully gotten past initialization and the JES initialization data is defined.
    def toJes: Option[JesBackendInitializationData] = genericInitializationData collectFirst { case d: JesBackendInitializationData => d }
  }

  val staticRuntimeAttributeDefinitions = {
    import JesRuntimeAttributes._
    import RuntimeAttributesKeys._

    Set(
      RuntimeAttributeDefinition(DockerKey, None, usedInCallCaching = true),
      RuntimeAttributeDefinition(ContinueOnReturnCodeKey, Option(staticDefaults(ContinueOnReturnCodeKey)), usedInCallCaching = true),
      RuntimeAttributeDefinition(CpuKey, Option(staticDefaults(CpuKey)), usedInCallCaching = false),
      RuntimeAttributeDefinition(FailOnStderrKey, Option(staticDefaults(FailOnStderrKey)), usedInCallCaching = true),
      RuntimeAttributeDefinition(MemoryKey, Option(staticDefaults(MemoryKey)), usedInCallCaching = false),
      RuntimeAttributeDefinition(DisksKey, Option(staticDefaults(DisksKey)), usedInCallCaching = false),
      RuntimeAttributeDefinition(ZonesKey, Option(staticDefaults(ZonesKey)), usedInCallCaching = false),
      RuntimeAttributeDefinition(PreemptibleKey, Option(staticDefaults(PreemptibleKey)), usedInCallCaching = false),
      RuntimeAttributeDefinition(BootDiskSizeKey, Option(staticDefaults(BootDiskSizeKey)), usedInCallCaching = false),
      RuntimeAttributeDefinition(NoAddressKey, Option(staticDefaults(NoAddressKey)), usedInCallCaching = false)
    )
  }
}
