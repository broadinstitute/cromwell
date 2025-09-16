package cromwell.services.healthmonitor.impl

import akka.actor.ActorRef
import com.typesafe.config.Config
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActor
import cromwell.services.healthmonitor.impl.common.EngineDatabaseMonitor

final class HealthMonitorServiceActor(val serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
    extends ProtoHealthMonitorServiceActor
    with EngineDatabaseMonitor {

  override lazy val subsystems: Set[ProtoHealthMonitorServiceActor.MonitoredSubsystem] = {

    val engineDatabaseSubsystemOption = if (serviceConfig.getBoolean("check-engine-database")) Some(EngineDb) else None

    Set(
      engineDatabaseSubsystemOption
    ).flatten
  }
}
