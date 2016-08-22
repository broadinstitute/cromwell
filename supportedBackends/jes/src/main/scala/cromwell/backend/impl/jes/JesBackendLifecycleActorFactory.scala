package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.{ActorRef, ActorSystem, Props}
import com.typesafe.config.Config
import cromwell.backend._
import cromwell.backend.callcaching.FileHasherWorkerActor.FileHashingFunction
import cromwell.backend.impl.jes.callcaching.JesBackendFileHashing
import cromwell.backend.validation.RuntimeAttributesKeys
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.{ExecutionStore, OutputStore}
import wdl4s.Call
import wdl4s.expression.WdlStandardLibraryFunctions

import scala.language.postfixOps

object JesBackendLifecycleActorFactory {
  implicit class Jessify(val genericInitializationData: Option[BackendInitializationData]) {
    // This leaves the result in an `Option` as finalization will be called even if initialization has failed, and if
    // initialization fails there won't be any initialization data.  The various `.get`s that occur below are in instances
    // where the workflow has successfully gotten past initialization and the JES initialization data is defined.
    def toJes: Option[JesBackendInitializationData] = genericInitializationData collectFirst { case d: JesBackendInitializationData => d }
  }
}

case class JesBackendLifecycleActorFactory(configurationDescriptor: BackendConfigurationDescriptor, actorSystem: ActorSystem) extends BackendLifecycleActorFactory {
  import JesBackendLifecycleActorFactory._

  val jesConfiguration = new JesConfiguration(configurationDescriptor)

  override def workflowInitializationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                                calls: Seq[Call],
                                                serviceRegistryActor: ActorRef): Option[Props] = {
    Option(JesInitializationActor.props(workflowDescriptor, calls, jesConfiguration, serviceRegistryActor).withDispatcher(BackendDispatcher))
  }

  override def jobExecutionActorProps(jobDescriptor: BackendJobDescriptor,
                                      initializationData: Option[BackendInitializationData],
                                      serviceRegistryActor: ActorRef): Props = {
    JesJobExecutionActor.props(jobDescriptor, jesConfiguration, initializationData.toJes.get, serviceRegistryActor).withDispatcher(BackendDispatcher)
  }

  override def workflowFinalizationActorProps(workflowDescriptor: BackendWorkflowDescriptor,
                                              calls: Seq[Call],
                                              executionStore: ExecutionStore,
                                              outputStore: OutputStore,
                                              initializationData: Option[BackendInitializationData]) = {
    Option(JesFinalizationActor.props(workflowDescriptor, calls, jesConfiguration, executionStore, outputStore, initializationData.toJes).withDispatcher(BackendDispatcher))
  }

  override val runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition] = {
    import RuntimeAttributesKeys._
    import JesRuntimeAttributes._

    Set(
      RuntimeAttributeDefinition(DockerKey, required = false, None, usedInCallCaching = true),
      RuntimeAttributeDefinition(ContinueOnReturnCodeKey, required = false, Some(staticDefaults(ContinueOnReturnCodeKey)), usedInCallCaching = true),
      RuntimeAttributeDefinition(CpuKey, required = false, Some(staticDefaults(CpuKey)), usedInCallCaching = false),
      RuntimeAttributeDefinition(FailOnStderrKey, required = false, Some(staticDefaults(FailOnStderrKey)), usedInCallCaching = true),
      RuntimeAttributeDefinition(MemoryKey, required = false, Some(staticDefaults(MemoryKey)), usedInCallCaching = false),
      RuntimeAttributeDefinition(DisksKey, required = false, Some(staticDefaults(DisksKey)), usedInCallCaching = false),
      RuntimeAttributeDefinition(ZonesKey, required = false, Some(staticDefaults(ZonesKey)), usedInCallCaching = false),
      RuntimeAttributeDefinition(PreemptibleKey, required = false, Some(staticDefaults(PreemptibleKey)), usedInCallCaching = false),
      RuntimeAttributeDefinition(BootDiskSizeKey, required = false, Some(staticDefaults(BootDiskSizeKey)), usedInCallCaching = false)
    )
  }

  override def expressionLanguageFunctions(workflowDescriptor: BackendWorkflowDescriptor,
                                           jobKey: BackendJobDescriptorKey,
                                           initializationData: Option[BackendInitializationData]): WdlStandardLibraryFunctions = {

    val jesCallPaths = initializationData.toJes.get.workflowPaths.toJesCallPaths(jobKey)
    new JesExpressionFunctions(List(jesCallPaths.gcsFileSystemWithUserAuth), jesCallPaths.callContext)
  }

  override def getExecutionRootPath(workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config,
                                    initializationData: Option[BackendInitializationData]): Path = {
    initializationData.toJes.get.workflowPaths.rootPath
  }

  override lazy val fileHashingFunction: Option[FileHashingFunction] = Option(FileHashingFunction(JesBackendFileHashing.getCrc32c))
}
