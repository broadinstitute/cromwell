package cromwell.backend.impl.sfs.config

import akka.actor.ActorSystem
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.callcaching.BackendHashingMethods
import cromwell.backend.sfs._
import lenthall.config.ScalaConfig._
import cromwell.backend.impl.sfs.config.ConfigConstants._
import cromwell.backend.sfs.callcaching.ConfigSfsBackendHashingMethods

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

  override val backendHashingMethods: BackendHashingMethods = ConfigSfsBackendHashingMethods(actorSystem)
}
