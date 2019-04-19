package cromwell.services.healthmonitor.impl.common

import cromwell.services.EngineServicesStore
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActor.{MonitoredSubsystem, SubsystemStatus, OkStatus}
import cats.instances.future._
import cats.syntax.functor._
import scala.concurrent.{ExecutionContext, Future}

/**
  * A mixin which provides a monitor to verify connectivity to the Cromwell engine database. Ideally this would be an
  * Object and not a Trait but the need for an ExecutionContext made that difficult with the () => Future[Subsystem]
  * signature of the checks.
  */
trait EngineDatabaseMonitor {
  implicit val ec: ExecutionContext

  lazy val EngineDb = MonitoredSubsystem("Engine Database", checkEngineDb _)

  /**
    * Demonstrates connectivity to the engine database by periodically making a small query
    */
  private def checkEngineDb(): Future[SubsystemStatus] = {
    EngineServicesStore.engineDatabaseInterface.queryDockerHashStoreEntries("DOESNOTEXIST") as OkStatus
  }
}
