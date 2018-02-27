package cromwell.services.healthmonitor

import akka.actor.Props
import akka.pattern.ask
import akka.util.Timeout
import com.typesafe.config.ConfigFactory
import cromwell.core.TestKitSuite
import cromwell.services.healthmonitor.HealthMonitorServiceActor.{GetCurrentStatus, MonitoredSubsystem, StatusCheckResponse, SubsystemStatus, _}
import cromwell.services.healthmonitor.HealthMonitorServiceActorSpec.{MockHealthMonitorServiceActor, _}
import org.scalatest.AsyncFlatSpecLike

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.language.postfixOps

class HealthMonitorServiceActorSpec extends TestKitSuite with AsyncFlatSpecLike {
  implicit val timeout = Timeout(5.seconds)

  "HealthMonitor" should "start with unknown status for all subsystems" in {
    newHealthMonitorActor().ask(GetCurrentStatus).mapTo[StatusCheckResponse] map { resp =>
      assert(!resp.ok)
      assert(resp.systems.filterNot(_._2.ok).filter(_._2 == UnknownStatus).keySet == MultipleSubsystems.map(_.name))
    }
  }

  it should "return ok status when all subsystems succeed" in {
    val hm = newHealthMonitorActor(Set(SuccessSubsystem))
    // Need to give a little bit of time for the cache to update
    Thread.sleep(1000)
    hm.ask(GetCurrentStatus).mapTo[StatusCheckResponse] map { resp =>
      assert(resp.ok)
      assert(resp.systems.collect({ case (name, status) if status.ok => name }).head == SuccessSubsystem.name)
    }
  }

  it should "fail if any subsystems fail but correctly bin them" in {
    val hm = newHealthMonitorActor()
    // Need to give a little bit of time for the cache to update
    Thread.sleep(1000)
    hm.ask(GetCurrentStatus).mapTo[StatusCheckResponse] map { resp =>
      assert(!resp.ok)
      val (ok, notOk) = resp.systems.partition { case (_, status) => status.ok }
      assert(ok.head._1 == SuccessSubsystem.name)
      assert(notOk.head._1 == FailSubsystem.name)
    }
  }

  it should "properly invalidate stale cache entries" in {
    val hm = newHealthMonitorActor()
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

    val timeoutSubsystem = MonitoredSubsystem("Timeout", timeOutCheck _)

    val hm = newHealthMonitorActor(Set(timeoutSubsystem))
    // Need to give a little bit of time for the cache to update
    Thread.sleep(1000)
    hm.ask(GetCurrentStatus).mapTo[StatusCheckResponse] map { resp =>
      assert(!resp.ok)
      assert(resp.systems.keySet.head == timeoutSubsystem.name)
      val errorMessages = resp.systems("Timeout").messages.get
      assert(errorMessages.size == 1)
      assert(errorMessages.head.startsWith("Timed out"))
    }
  }

  it should "record the error messages of the failures" in {
    val hm = newHealthMonitorActor(Set(FailSubsystem))
    // Need to give a little bit of time for the cache to update
    Thread.sleep(1000)
    hm.ask(GetCurrentStatus).mapTo[StatusCheckResponse] map { resp =>
      assert(!resp.ok)
      assert(resp.systems.keySet.head == FailSubsystem.name)
      val errorMessages = resp.systems("Failure").messages.get
      assert(errorMessages.size == 1)
      assert(errorMessages.head == "womp womp")
    }
  }

  it should "retry failures the appropriate amount of times before failing" in {
    // Start up an actor with a status check function that's ok the first time but subsequently fails, incrementing a
    // counter each time it is invoked.
    var statusCheckCount = 0
    def okThenFailStatusChecker(): Future[SubsystemStatus] = {
      statusCheckCount = statusCheckCount + 1
      if (statusCheckCount == 1) {
        Future.successful(SubsystemStatus(ok = true, messages = None))
      } else {
        Future.failed(new RuntimeException("womp womp"))
      }
    }

    var storeInvocationCount = 0
    // Check the status that's being recorded, incrementing a counter each time it is invoked.
    var statusStores = List.empty[SubsystemStatus]
    def storeChecker(monitoredSubsystem: MonitoredSubsystem, subsystemStatus: SubsystemStatus): Unit = {
      storeInvocationCount = storeInvocationCount + 1
      assert(subsystemStatus.ok ^ (statusCheckCount == 4))
      statusStores = statusStores :+ subsystemStatus
    }

    val failAndCountSubsystem = MonitoredSubsystem("FailAndCountSubsystem", okThenFailStatusChecker _)
    val hm = newHealthMonitorActor(Set(failAndCountSubsystem), retryCount = 2, storeChecker = Option(storeChecker))
    hm ! Check(failAndCountSubsystem, 2)

    Thread.sleep(2000L)
    // Assert status was stored exactly twice, initially ok and then not ok.
    assert(statusStores.map(_.ok) == List(true, false))
    Future.successful(assert(statusCheckCount == 4))
  }

  private def newHealthMonitorActor(subsystems: Set[MonitoredSubsystem] = MultipleSubsystems, retryCount: Int = 0, storeChecker: Option[(MonitoredSubsystem, SubsystemStatus) => Unit] = None) = {
    system.actorOf(Props(new MockHealthMonitorServiceActor(subsystems, retryCount, storeChecker)))
  }
}

object HealthMonitorServiceActorSpec {
  val SuccessSubsystem = MonitoredSubsystem("Success", mockCheckSuccess _)
  val FailSubsystem = MonitoredSubsystem("Failure", mockCheckFailure _)
  val MultipleSubsystems = Set(SuccessSubsystem, FailSubsystem)

  class MockHealthMonitorServiceActor(override val subsystems: Set[MonitoredSubsystem], retryCount: Int = 0, storeChecker: Option[(MonitoredSubsystem, SubsystemStatus) => Unit] = None) extends HealthMonitorServiceActor {
    override val serviceConfig = ConfigFactory.empty
    override lazy val futureTimeout = 500 milliseconds
    override lazy val staleThreshold = 3 seconds
    override lazy val sweepInterval = 1 second
    override lazy val failureRetryInterval = 200 milliseconds
    override lazy val failureRetryCount = retryCount

    override private[healthmonitor] def scheduleSweepCheck(subsystem: MonitoredSubsystem): Unit = ()

    def statusFor(monitoredSubsystem: MonitoredSubsystem): CachedSubsystemStatus = statusCache(monitoredSubsystem)

    override def store(monitoredSubsystem: MonitoredSubsystem, subsystemStatus: SubsystemStatus): Unit = {
      storeChecker foreach (_.apply(monitoredSubsystem, subsystemStatus))
      super.store(monitoredSubsystem, subsystemStatus)
    }
  }

  private def mockCheckSuccess(): Future[SubsystemStatus] = {
    Future.successful(OkStatus)
  }

  private def mockCheckFailure(): Future[SubsystemStatus] = {
    Future.failed(new RuntimeException("womp womp"))
  }
}