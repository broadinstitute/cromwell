package cromwell.services.healthmonitor

import akka.actor.ActorRef
import com.typesafe.config.Config
import cromwell.services.healthmonitor.impl.common.{DockerHubMonitor, EngineDatabaseMonitor}

final class ChooseYourOwnAdventureHealthMonitorServiceActor() (val serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef)
  extends HealthMonitorServiceActor
    with DockerHubMonitor
    with EngineDatabaseMonitor {
  override implicit val system = context.system





  override def subsystems: Set[HealthMonitorServiceActor.MonitoredSubsystem] = ???
}
