package cromwell.services.healthmonitor

import akka.actor.{ActorRef, Props}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.TestKitSuite
import cromwell.services.healthmonitor.HealthMonitorServiceActor.{GetCurrentStatus, MonitoredSubsystem, StatusCheckResponse, SubsystemStatus}
import akka.pattern.ask

import scala.concurrent.duration._
import HealthMonitorServiceActor._
import HealthMonitorServiceActorSpec._
import akka.util.Timeout
import cromwell.services.healthmonitor.HealthMonitorServiceActorSpec.MockHealthMonitorServiceActor
import org.scalatest.AsyncFlatSpecLike

import scala.concurrent.Future

class HealthMonitorServiceActorSpec extends TestKitSuite with AsyncFlatSpecLike {
  implicit val timeout = Timeout(5.seconds)

  "HealthMonitor" should "start with unknown status for all subsystems" in {
    newHealthMonitorActor().ask(GetCurrentStatus).mapTo[StatusCheckResponse] map { resp =>
      assert(!resp.ok)
      assert(resp.systems.filterNot(_._2.ok).filter(_._2 == UnknownStatus).keySet == MultipleSubsystems.map(_.name))
    }
  }

  it should "return ok status when all subsystems succeed" in {
    val hm = newHealthMonitorActor(SuccessSubsystem)
    hm ! CheckAll
    // Need to give a little bit of time for the cache to update
    Thread.sleep(1000)
    hm.ask(GetCurrentStatus).mapTo[StatusCheckResponse] map { resp =>
      assert(resp.ok)
      assert(resp.systems.filter(_._2.ok).keySet == SuccessSubsystem.map(_.name))
    }
  }

  it should "fail if any subsystems fail but correctly bin them" in {
    val hm = newHealthMonitorActor()
    hm ! CheckAll
    // Need to give a little bit of time for the cache to update
    Thread.sleep(1000)
    hm.ask(GetCurrentStatus).mapTo[StatusCheckResponse] map { resp =>
      assert(!resp.ok)
      assert(resp.systems.filter(_._2.ok).keySet == SuccessSubsystem.map(_.name))
      assert(resp.systems.filterNot(_._2.ok).filterNot(_._2 == UnknownStatus).keySet == FailSubsystem.map(_.name))
    }
  }

  it should "properly invalidate stale cache entries" in {
    val hm = newHealthMonitorActor()
    hm ! CheckAll
    // Need to give a little bit of time for the cache to update and then go stale
    Thread.sleep(5000)
    hm.ask(GetCurrentStatus).mapTo[StatusCheckResponse] map { resp =>
      assert(!resp.ok)
      assert(resp.systems.filterNot(_._2.ok).filter(_._2 == UnknownStatus).keySet == MultipleSubsystems.map(_.name))
    }

  }

  it should "handle timed out futures" in {
    def timeOutCheck(): Future[SubsystemStatus] = {
      Future {
        Thread.sleep(3000)
        OkStatus
      }
    }

    val timeoutSubsystem = Set(MonitoredSubsystem("Timeout", timeOutCheck _))

    val hm = newHealthMonitorActor(timeoutSubsystem)
    hm ! CheckAll
    // Need to give a little bit of time for the cache to update
    Thread.sleep(1000)
    hm.ask(GetCurrentStatus).mapTo[StatusCheckResponse] map { resp =>
      assert(!resp.ok)
      assert(resp.systems.filterNot(_._2.ok).filterNot(_._2 == UnknownStatus).keySet == timeoutSubsystem.map(_.name))
      val errorMessages = resp.systems("Timeout").messages.get
      assert(errorMessages.size == 1)
      assert(errorMessages(0).startsWith("Timed out"))
    }
  }

  it should "record the error messages of the failures" in {
    val hm = newHealthMonitorActor(FailSubsystem)
    hm ! CheckAll
    // Need to give a little bit of time for the cache to update
    Thread.sleep(1000)
    hm.ask(GetCurrentStatus).mapTo[StatusCheckResponse] map { resp =>
      assert(!resp.ok)
      assert(resp.systems.filterNot(_._2.ok).filterNot(_._2 == UnknownStatus).keySet == FailSubsystem.map(_.name))
      val errorMessages = resp.systems("Failure").messages.get
      assert(errorMessages.size == 1)
      assert(errorMessages(0) == "womp womp")

    }
  }

  private def newHealthMonitorActor(subsystems: Set[MonitoredSubsystem] = MultipleSubsystems): ActorRef = {
    system.actorOf(Props(new MockHealthMonitorServiceActor(ConfigFactory.empty, subsystems)))
  }
}

object HealthMonitorServiceActorSpec {
  val SuccessSubsystem = Set(MonitoredSubsystem("Success", mockCheckSuccess _))
  val FailSubsystem = Set(MonitoredSubsystem("Failure", mockCheckFailure _))
  val MultipleSubsystems = SuccessSubsystem ++ FailSubsystem


  class MockHealthMonitorServiceActor(val serviceConfig: Config, override val subsystems: Set[MonitoredSubsystem]) extends HealthMonitorServiceActor {
    override lazy val futureTimeout = 500.milli
    override lazy val staleThreshold = 3.seconds
  }

  private def mockCheckSuccess(): Future[SubsystemStatus] = {
    Future.successful(OkStatus)
  }

  private def mockCheckFailure(): Future[SubsystemStatus] = {
    Future.failed(new RuntimeException("womp womp"))
  }
}