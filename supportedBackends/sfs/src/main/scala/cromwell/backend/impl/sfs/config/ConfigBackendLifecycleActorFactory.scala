package cromwell.backend.impl.sfs.config

import com.typesafe.config.Config
import cromwell.backend.callcaching.FileHashingActor.FileHashingFunction
import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.sfs._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, RuntimeAttributeDefinition}
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import net.ceedubs.ficus.Ficus._
import org.slf4j.LoggerFactory

/**
  * Builds a backend by reading the job control from the config.
  *
  * @param configurationDescriptor The config information.
  */
class ConfigBackendLifecycleActorFactory(name: String, val configurationDescriptor: BackendConfigurationDescriptor)
  extends SharedFileSystemBackendLifecycleActorFactory {

  lazy val logger = LoggerFactory.getLogger(getClass)
  lazy val hashingStrategy = {
    configurationDescriptor.backendConfig.as[Option[Config]]("filesystems.local.caching") map ConfigHashingStrategy.apply getOrElse ConfigHashingStrategy.defaultStrategy
  }

  override def initializationActorClass = classOf[ConfigInitializationActor]

  override def asyncJobExecutionActorClass: Class[_ <: ConfigAsyncJobExecutionActor] = {
    val runInBackground = configurationDescriptor.backendConfig.as[Option[Boolean]](RunInBackgroundConfig).getOrElse(false)
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

  override lazy val fileHashingFunction: Option[FileHashingFunction] = {
    logger.debug(hashingStrategy.toString)
    Option(FileHashingFunction(hashingStrategy.getHash))
  }

  override lazy val fileHashingActorCount: Int = 5

  override val jobExecutionTokenType: JobExecutionTokenType = {
    val concurrentJobLimit = configurationDescriptor.backendConfig.as[Option[Int]]("concurrent-job-limit")
    JobExecutionTokenType(name, concurrentJobLimit)
  }
}
