package cromwell.services.healthmonitor

import java.util.concurrent.TimeoutException

import akka.actor.{Actor, Scheduler, Timers}
import akka.pattern.{after, pipe}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import com.typesafe.scalalogging.LazyLogging

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._
import scala.util.control.NonFatal
import scala.language.postfixOps
import HealthMonitorServiceActor._
import com.typesafe.config.Config
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import net.ceedubs.ficus.Ficus._

/**
  * A trait which provides a health monitoring service. Implementers are expected to provide a collection of functions
  * which will each check some underlying service on which Cromwell depends, for instance connectivity to the database
  * or DockerHub. At a regular interval the health monitor actor will run each of these functions and keep their
  * current state in a cache.
  *
  * This actor may also be interrogated for the current status of the monitored subsystems. The response will indicate
  * if all of the subsystems are deemed to be healthy, and if this is not true subsystem-specific messages will be
  * returned to the client.
  */
trait HealthMonitorServiceActor extends Actor with LazyLogging with Timers {
  val serviceConfig: Config
  def subsystems: Set[MonitoredSubsystem]

  implicit val ec: ExecutionContext = context.system.dispatcher

  lazy val futureTimeout: FiniteDuration = serviceConfig.as[Option[FiniteDuration]]("services.HealthMonitor.check-timeout").getOrElse(DefaultFutureTimeout)
  lazy val staleThreshold: FiniteDuration = serviceConfig.as[Option[FiniteDuration]]("services.HealthMonitor.status-ttl").getOrElse(DefaultStaleThreshold)

  override def preStart(): Unit = {
    val sweepTime = serviceConfig.as[Option[FiniteDuration]]("services.HealthMonitor.check-refresh-time").getOrElse(DefaultSweepTime)
    timers.startPeriodicTimer(CheckTickKey, CheckAll, sweepTime)
    logger.info("Starting health monitor with the following checks: " + subsystems.map(_.name).mkString(", "))
  }

  /**
    * Contains each subsystem status along with a timestamp of when the entry was made so we know when the status
    * goes stale. Initialized with unknown status.
    */
  private var statusCache: Map[MonitoredSubsystem, CachedSubsystemStatus] = {
    val now = System.currentTimeMillis
    subsystems.map((_, CachedSubsystemStatus(UnknownStatus, now))).toMap
  }

  override def receive: Receive = {
    case CheckAll =>
      logger.debug("Checking status on all subsystems")
      subsystems.foreach(checkSubsystem)
    case Store(subsystem, status) => store(subsystem, status)
    case GetCurrentStatus => sender ! getCurrentStatus
    case ShutdownCommand => context.stop(self) // Not necessary but service registry requires it. See #2575
  }

  private def checkSubsystem(subsystem: MonitoredSubsystem): Unit = {
    val result = subsystem.check()
    val errMsg = s"Timed out after ${futureTimeout.toString} waiting for a response from ${subsystem.toString}"
    result.withTimeout(futureTimeout, errMsg, context.system.scheduler)
      .recover { case NonFatal(ex) =>
        failedStatus(ex.getMessage)
      } map {
      Store(subsystem, _)
    } pipeTo self

    ()
  }

  private def store(subsystem: MonitoredSubsystem, status: SubsystemStatus): Unit = {
    statusCache = statusCache + (subsystem -> CachedSubsystemStatus(status, System.currentTimeMillis))
    logger.debug(s"New health monitor state: $statusCache")
  }

  private def getCurrentStatus: StatusCheckResponse = {
    val now = System.currentTimeMillis()
    // Convert any expired statuses to unknown
    val processed = statusCache map {
      case (s, c) if now - c.created > staleThreshold.toMillis => s.name -> UnknownStatus
      case (s, c) => s.name -> c.status
    }

    // overall status is ok iff all subsystems are ok
    val overall = processed.forall(_._2.ok)

    StatusCheckResponse(overall, processed)
  }
}

object HealthMonitorServiceActor {
  val DefaultFutureTimeout: FiniteDuration = 1 minute
  val DefaultStaleThreshold: FiniteDuration = 15 minutes
  val DefaultSweepTime: FiniteDuration = 1 minute // FIXME 5 minutes

  private case object CheckTickKey

  val OkStatus = SubsystemStatus(true, None)
  val UnknownStatus = SubsystemStatus(false, Option(List("Unknown status")))
  def failedStatus(message: String) = SubsystemStatus(false, Option(List(message)))

  final case class MonitoredSubsystem(name: String, check: () => Future[SubsystemStatus])
  final case class SubsystemStatus(ok: Boolean, messages: Option[List[String]])
  final case class CachedSubsystemStatus(status: SubsystemStatus, created: Long) // created is time in millis when status was captured

  sealed abstract class HealthMonitorServiceActorRequest
  case object CheckAll extends HealthMonitorServiceActorRequest
  final case class Store(subsystem: MonitoredSubsystem, status: SubsystemStatus) extends HealthMonitorServiceActorRequest
  case object GetCurrentStatus extends HealthMonitorServiceActorRequest with ServiceRegistryMessage { override val serviceName = "HealthMonitor" }

  sealed abstract class HealthMonitorServiceActorResponse
  final case class StatusCheckResponse(ok: Boolean, systems: Map[String, SubsystemStatus]) extends HealthMonitorServiceActorResponse

  /**
    * Adds non-blocking timeout support to futures.
    * Example usage:
    * {{{
    *   val future = Future(Thread.sleep(1000*60*60*24*365)) // 1 year
    *   Await.result(future.withTimeout(5 seconds, "Timed out"), 365 days)
    *   // returns in 5 seconds
    * }}}
    */
  private implicit class FutureWithTimeout[A](f: Future[A]) {
    def withTimeout(duration: FiniteDuration, errMsg: String, scheduler: Scheduler)(implicit ec: ExecutionContext): Future[A] =
      Future.firstCompletedOf(List(f, after(duration, scheduler)(Future.failed(new TimeoutException(errMsg)))))
  }
}
