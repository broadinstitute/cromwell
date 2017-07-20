package cromwell.services.healthmonitor.impl.standard

import com.typesafe.config.Config
import cromwell.services.healthmonitor.HealthMonitorServiceActor
import cromwell.services.healthmonitor.HealthMonitorServiceActor.MonitoredSubsystem
import cromwell.services.healthmonitor.impl.CommonMonitoredSubsystems

/**
  * A health monitor implementation which only monitors things available to out of the box Cromwell installations
  */
class StandardHealthMonitorServiceActor(val serviceConfig: Config, globalConfig: Config)
  extends HealthMonitorServiceActor
    with CommonMonitoredSubsystems {
  override lazy val subsystems: List[MonitoredSubsystem] = List(DockerHub, EngineDb)
}
