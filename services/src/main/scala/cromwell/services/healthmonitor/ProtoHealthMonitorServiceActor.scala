package cromwell.services.healthmonitor

import java.util.UUID
import java.util.concurrent.TimeoutException

import akka.actor.{Actor, Scheduler, Status, Timers}
import akka.pattern.after
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActor._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util._
import scala.util.control.NonFatal

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
trait ProtoHealthMonitorServiceActor extends Actor with LazyLogging with Timers {
  val serviceConfig: Config
  def subsystems: Set[MonitoredSubsystem]

  implicit val ec: ExecutionContext = context.system.dispatcher

  lazy val sweepInterval = serviceConfig.as[Option[FiniteDuration]]("services.HealthMonitor.check-refresh-time").getOrElse(DefaultSweepTime)
  lazy val futureTimeout: FiniteDuration = serviceConfig.as[Option[FiniteDuration]]("services.HealthMonitor.check-timeout").getOrElse(DefaultFutureTimeout)
  lazy val staleThreshold: FiniteDuration = serviceConfig.as[Option[FiniteDuration]]("services.HealthMonitor.status-ttl").getOrElse(DefaultStaleThreshold)
  lazy val failureRetryCount: Int = serviceConfig.as[Option[Int]]("services.HealthMonitor.check-failure-retry-count").getOrElse(DefaultFailureRetryCount)
  lazy val failureRetryInterval: FiniteDuration = serviceConfig.as[Option[FiniteDuration]]("services.HealthMonitor.check-failure-retry-interval").getOrElse(DefaultFailureRetryInterval)

  /**
    * Contains each subsystem status along with a timestamp of when the entry was made so we know when the status
    * goes stale. Initialized with unknown status.
    */
  private[healthmonitor] var statusCache: Map[MonitoredSubsystem, CachedSubsystemStatus] = {
    val now = System.currentTimeMillis
    subsystems.map((_, CachedSubsystemStatus(UnknownStatus, now))).toMap
  }

  private[healthmonitor] def initialize(): Unit = {
    subsystems foreach { s =>
      logger.info(s"Availability of '${s.name}' will be monitored and reported via the '/engine/v1/status' API")
      self ! Check(s, failureRetryCount)
    }
  }

  private def check(subsystem: MonitoredSubsystem, after: FiniteDuration, withRetriesLeft: Int): Unit = {
    timers.startSingleTimer(UUID.randomUUID(), Check(subsystem, withRetriesLeft), after)
  }

  private[healthmonitor] def scheduleFailedRetryCheck(subsystem: MonitoredSubsystem, retriesLeft: Int): Unit = {
    check(subsystem, after = failureRetryInterval, withRetriesLeft = retriesLeft - 1)
  }

  private[healthmonitor] def scheduleSweepCheck(subsystem: MonitoredSubsystem): Unit = {
    check(subsystem, after = sweepInterval, withRetriesLeft = failureRetryCount)
  }

  initialize()

  override def receive: Receive = {
    case Check(subsystem, retriesLeft) =>
      checkSubsystem(subsystem, retriesLeft)
    case Store(subsystem, status) =>
      store(subsystem, status)
      scheduleSweepCheck(subsystem)
    case GetCurrentStatus => sender ! getCurrentStatus
    case ShutdownCommand => context.stop(self) // Not necessary but service registry requires it. See #2575
    case Status.Failure(f) => logger.error("Unexpected Status.Failure received", f)
    case e => logger.error("Unexpected Status.Failure received: {}", e.toString)
  }

  private def checkSubsystem(subsystem: MonitoredSubsystem, retriesLeft: Int): Unit = {
    val result = subsystem.check()
    val errMsg = s"Timed out after ${futureTimeout.toString} waiting for a response from ${subsystem.toString}"
    result.withTimeout(futureTimeout, errMsg, context.system.scheduler) onComplete {
      case Success(status) if status.ok =>
        self ! Store(subsystem, status)
      case Success(_) if retriesLeft > 0 =>
        // This is a success in getting the status but the status is a failure. However there are retries available.
        scheduleFailedRetryCheck(subsystem, retriesLeft)
      case Success(status) =>
        // This is also a success in getting the status but the status is a failure and there are no more retries available.
        // Womp womp.
        self ! Store(subsystem, status)
      case Failure(NonFatal(_)) if retriesLeft > 0 =>
        scheduleFailedRetryCheck(subsystem, retriesLeft)
      case Failure(NonFatal(e)) =>
        self ! Store(subsystem, SubsystemStatus(ok = false, messages = Option(List(e.getMessage))))
      case f =>
        self ! f
    }
    ()
  }

  private[healthmonitor] def store(subsystem: MonitoredSubsystem, status: SubsystemStatus): Unit = {
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

object ProtoHealthMonitorServiceActor {
  val DefaultFutureTimeout: FiniteDuration = 1 minute
  val DefaultStaleThreshold: FiniteDuration = 15 minutes
  val DefaultSweepTime: FiniteDuration = 5 minutes
  val DefaultFailureRetryCount: Int = 3
  val DefaultFailureRetryInterval: FiniteDuration = 30 seconds

  val OkStatus = SubsystemStatus(true, None)
  val UnknownStatus = SubsystemStatus(false, Option(List("Unknown status")))
  def failedStatus(message: String) = SubsystemStatus(false, Option(List(message)))

  final case class MonitoredSubsystem(name: String, check: () => Future[SubsystemStatus])
  final case class SubsystemStatus(ok: Boolean, messages: Option[List[String]])
  final case class CachedSubsystemStatus(status: SubsystemStatus, created: Long) // created is time in millis when status was captured

  sealed abstract class HealthMonitorServiceActorRequest
  case class Check(subsystem: MonitoredSubsystem, retriesLeft: Int) extends HealthMonitorServiceActorRequest
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
