package cromwell.backend.impl.sfs.config

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.sfs._
import lenthall.config.ScalaConfig._
import cromwell.backend.impl.sfs.config.ConfigConstants._

/**
  * Builds a backend by reading the job control from the config.
  *
  * @param configurationDescriptor The config information.
  */
class ConfigBackendLifecycleActorFactory(override val configurationDescriptor: BackendConfigurationDescriptor)
  extends SharedFileSystemBackendLifecycleActorFactory {

  override def initializationActorClass = classOf[ConfigInitializationActor]

  override def asyncJobExecutionActorClass: Class[_ <: ConfigAsyncJobExecutionActor] = {
    val runInBackground = configurationDescriptor.backendConfig.getBooleanOr(RunInBackgroundConfig, default = false)
    if (runInBackground)
      classOf[BackgroundConfigAsyncJobExecutionActor]
    else
      classOf[DispatchedConfigAsyncJobExecutionActor]
  }
}
