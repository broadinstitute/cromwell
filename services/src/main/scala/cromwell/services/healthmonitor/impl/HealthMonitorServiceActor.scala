package cromwell.services.healthmonitor.impl

import akka.actor.ActorRef
import com.typesafe.config.Config
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActor
import cromwell.services.healthmonitor.impl.workbench.WorkbenchHealthMonitorServiceActor

final class HealthMonitorServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
  extends WorkbenchHealthMonitorServiceActor(serviceConfig, globalConfig, serviceRegistryActor) {
  override implicit val system = context.system

  override lazy val subsystems: Set[ProtoHealthMonitorServiceActor.MonitoredSubsystem] = {

    val dockerHubSubsystemOption = if (serviceConfig.getBoolean("check-dockerhub")) Some(DockerHub) else None
    val engineDatabaseSubsystemOption = if (serviceConfig.getBoolean("check-engine-database")) Some(EngineDb) else None
    val gcsSubsystemOption = if (serviceConfig.getBoolean("check-gcs")) Some(Gcs) else None

    Set(
      dockerHubSubsystemOption,
      engineDatabaseSubsystemOption,
      gcsSubsystemOption,
    ).flatten ++ PapiSubsystems
  }
}
