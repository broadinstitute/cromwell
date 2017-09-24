package cromwell.services.healthmonitor.impl.standard

import com.typesafe.config.Config
import cromwell.services.healthmonitor.HealthMonitorServiceActor
import cromwell.services.healthmonitor.HealthMonitorServiceActor.MonitoredSubsystem
import cromwell.services.healthmonitor.impl.common.{DockerHubMonitor, EngineDatabaseMonitor}

/**
  * A health monitor implementation which only monitors things available to out of the box Cromwell installations
  */
class StandardHealthMonitorServiceActor(val serviceConfig: Config, globalConfig: Config)
  extends HealthMonitorServiceActor
    with DockerHubMonitor
    with EngineDatabaseMonitor {
  override implicit val system = context.system
  override lazy val subsystems: Set[MonitoredSubsystem] = Set(DockerHub, EngineDb)
}
