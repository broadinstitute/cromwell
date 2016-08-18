package cromwell.backend.impl.local

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import cromwell.backend.BackendConfigurationDescriptor
import cromwell.backend.impl.local.LocalBackendLifecycleActorFactory._
import cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory

@deprecated("Remains until travis/centaur is updated to stop using this class.", "SFS")
class LocalBackendLifecycleActorFactory(configurationDescriptor: BackendConfigurationDescriptor, actorSystem: ActorSystem)
  extends ConfigBackendLifecycleActorFactory(reconfig(configurationDescriptor), actorSystem)


@deprecated("Remains until travis/centaur is updated to stop using this class.", "SFS")
object LocalBackendLifecycleActorFactory {
  def reconfig(configurationDescriptor: BackendConfigurationDescriptor): BackendConfigurationDescriptor = {
    val backendConfig = configurationDescriptor.backendConfig
    val globalConfig = configurationDescriptor.globalConfig
    BackendConfigurationDescriptor(backendConfig.withFallback(localConfig), globalConfig)
  }

  val localConfig = ConfigFactory.parseString(
    """
      |run-in-background = true
      |runtime-attributes = "String? docker"
      |submit = "/bin/bash ${script}"
      |submit-docker = "docker run --rm -v ${cwd}:${docker_cwd} -i ${docker} /bin/bash < ${script}"
    """.stripMargin)
}
