package cromwell.services.healthmonitor.impl.noop

import com.typesafe.config.Config
import cromwell.services.healthmonitor.HealthMonitorServiceActor
import cromwell.services.healthmonitor.HealthMonitorServiceActor.MonitoredSubsystem

/**
  * A health monitor which performs no checks
  */
class NoopHealthMonitorServiceActor (val serviceConfig: Config, globalConfig: Config)
  extends HealthMonitorServiceActor {
  override val subsystems: Set[MonitoredSubsystem] = Set.empty[MonitoredSubsystem]
}
