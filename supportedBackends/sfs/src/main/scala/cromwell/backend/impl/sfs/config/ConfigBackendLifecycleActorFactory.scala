package cromwell.backend.impl.sfs.config

import cromwell.backend.callcaching.FileHashingActor.FileHashingFunction
import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.sfs._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, RuntimeAttributeDefinition}
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import lenthall.config.ScalaConfig._

/**
  * Builds a backend by reading the job control from the config.
  *
  * @param configurationDescriptor The config information.
  */
class ConfigBackendLifecycleActorFactory(name: String, val configurationDescriptor: BackendConfigurationDescriptor)
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
  override lazy val fileHashingActorCount: Int = 5

  override val jobExecutionTokenType: JobExecutionTokenType = {
    val concurrentJobLimit = configurationDescriptor.backendConfig.getIntOption("concurrent-job-limit")
    System.out.println(s"Operating with a concurrent job limit for $name of ${concurrentJobLimit.getOrElse("None")}")
    JobExecutionTokenType(name, concurrentJobLimit)
  }
}
