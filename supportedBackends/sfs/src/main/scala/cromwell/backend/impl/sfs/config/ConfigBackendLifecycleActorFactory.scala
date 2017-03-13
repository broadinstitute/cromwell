package cromwell.backend.impl.sfs.config

import com.typesafe.config.Config
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.sfs._
import cromwell.backend.standard.callcaching.StandardFileHashingActor
import cromwell.core.JobExecutionToken.JobExecutionTokenType
import net.ceedubs.ficus.Ficus._
import org.slf4j.{Logger, LoggerFactory}

/**
  * Builds a backend by reading the job control from the config.
  *
  * @param configurationDescriptor The config information.
  */
class ConfigBackendLifecycleActorFactory(name: String, val configurationDescriptor: BackendConfigurationDescriptor)
  extends SharedFileSystemBackendLifecycleActorFactory {

  lazy val logger: Logger = LoggerFactory.getLogger(getClass)
  lazy val hashingStrategy: ConfigHashingStrategy = {
    configurationDescriptor.backendConfig.as[Option[Config]]("filesystems.local.caching") map ConfigHashingStrategy.apply getOrElse ConfigHashingStrategy.defaultStrategy
  }

  override lazy val initializationActorClass: Class[ConfigInitializationActor] = classOf[ConfigInitializationActor]

  override lazy val asyncExecutionActorClass: Class[_ <: ConfigAsyncJobExecutionActor] = {
    val runInBackground = configurationDescriptor.backendConfig.as[Option[Boolean]](RunInBackgroundConfig).getOrElse(false)
    if (runInBackground)
      classOf[BackgroundConfigAsyncJobExecutionActor]
    else
      classOf[DispatchedConfigAsyncJobExecutionActor]
  }

  override lazy val fileHashingActorClassOption: Option[Class[_ <: StandardFileHashingActor]] = Option(classOf[ConfigBackendFileHashingActor])

  override val jobExecutionTokenType: JobExecutionTokenType = {
    val concurrentJobLimit = configurationDescriptor.backendConfig.as[Option[Int]]("concurrent-job-limit")
    JobExecutionTokenType(name, concurrentJobLimit)
  }
}
