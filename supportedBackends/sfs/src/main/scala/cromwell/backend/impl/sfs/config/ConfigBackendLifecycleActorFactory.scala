package cromwell.backend.impl.sfs.config

import cromwell.backend.callcaching.FileHasherWorkerActor.FileHashingFunction
import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.sfs._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, RuntimeAttributeDefinition}
import lenthall.config.ScalaConfig._

/**
  * Builds a backend by reading the job control from the config.
  *
  * @param configurationDescriptor The config information.
  */
class ConfigBackendLifecycleActorFactory(val configurationDescriptor: BackendConfigurationDescriptor)
  extends SharedFileSystemBackendLifecycleActorFactory {

  override def initializationActorClass = classOf[ConfigInitializationActor]

  override def asyncJobExecutionActorClass: Class[_ <: ConfigAsyncJobExecutionActor] = {
    val runInBackground = configurationDescriptor.backendConfig.getBooleanOr(RunInBackgroundConfig, default = false)
    if (runInBackground)
      classOf[BackgroundConfigAsyncJobExecutionActor]
    else
      classOf[DispatchedConfigAsyncJobExecutionActor]
  }

  override def runtimeAttributeDefinitions(initializationDataOption: Option[BackendInitializationData]):
  Set[RuntimeAttributeDefinition] = {
    val initializationData = BackendInitializationData.
      as[SharedFileSystemBackendInitializationData](initializationDataOption)

    initializationData.runtimeAttributesBuilder.definitions.toSet
  }

  override lazy val fileHashingFunction: Option[FileHashingFunction] = Option(FileHashingFunction(ConfigBackendFileHashing.getMd5Result))
  override lazy val fileHashingWorkerCount: Int = 5
}
