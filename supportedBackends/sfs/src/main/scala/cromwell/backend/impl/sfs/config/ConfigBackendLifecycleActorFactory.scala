package cromwell.backend.impl.sfs.config

import akka.actor.ActorSystem
import cromwell.backend.callcaching.FileHasherWorkerActor.FileHashingFunction
import cromwell.backend.{BackendConfigurationDescriptor, RuntimeAttributeDefinition}
import cromwell.backend.sfs._
import lenthall.config.ScalaConfig._
import cromwell.backend.impl.sfs.config.ConfigConstants._
import wdl4s.values.WdlBoolean

/**
  * Builds a backend by reading the job control from the config.
  *
  * @param configurationDescriptor The config information.
  */
class ConfigBackendLifecycleActorFactory(val configurationDescriptor: BackendConfigurationDescriptor, val actorSystem: ActorSystem)
  extends SharedFileSystemBackendLifecycleActorFactory {

  override def initializationActorClass = classOf[ConfigInitializationActor]

  override def asyncJobExecutionActorClass: Class[_ <: ConfigAsyncJobExecutionActor] = {
    val runInBackground = configurationDescriptor.backendConfig.getBooleanOr(RunInBackgroundConfig, default = false)
    if (runInBackground)
      classOf[BackgroundConfigAsyncJobExecutionActor]
    else
      classOf[DispatchedConfigAsyncJobExecutionActor]
  }

  // TODO: Wire in the funky configuration options too
  override def runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition] = {
    import cromwell.backend.validation.RuntimeAttributesKeys._
    Set(
      RuntimeAttributeDefinition(DockerKey, required = false, None, usedInCallCaching = true),
      RuntimeAttributeDefinition(ContinueOnReturnCodeKey, required = false, Some(WdlBoolean(false)), usedInCallCaching = true),
      RuntimeAttributeDefinition(CpuKey, required = false, None, usedInCallCaching = false),
      RuntimeAttributeDefinition(FailOnStderrKey, required = false, Some(WdlBoolean(true)), usedInCallCaching = true),
      RuntimeAttributeDefinition(MemoryKey, required = false, None, usedInCallCaching = false)
    )
  }

  override lazy val fileHashingFunction: Option[FileHashingFunction] = Option(FileHashingFunction(ConfigBackendFileHashing.getMd5Result))
  override lazy val fileHashingWorkerCount: Int = 5
}
