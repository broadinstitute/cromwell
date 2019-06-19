package cromwell.services.healthmonitor

import java.util.concurrent.Executors

import akka.actor.ActorRef
import akka.pattern.ask
import akka.testkit.TestActorRef
import akka.util.Timeout
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.TestKitSuite
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActor._
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActorSpec._
import org.scalatest.{Assertion, FlatSpecLike}
import org.scalatest.concurrent.{Eventually, ScaledTimeSpans}

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, ExecutionContextExecutor, Future}
import scala.language.postfixOps

class ProtoHealthMonitorServiceActorSpec extends TestKitSuite with FlatSpecLike with Eventually {
  implicit val timeout = Timeout(scaled(5.seconds))
  final implicit val blockingEc: ExecutionContextExecutor = ExecutionContext.fromExecutor(
    Executors.newCachedThreadPool()
  )

  override val patienceConfig = PatienceConfig(timeout = scaled(20 seconds), interval = scaled(1 second))
  implicit val patience = patienceConfig

  private def eventualStatus(actorRef: ActorRef, ok: Boolean, systems: (MonitoredSubsystem, SubsystemStatus)*): Assertion = {
    case class NameAndStatus(name: String, status: SubsystemStatus)

    eventually {
      val resp = Await.result(actorRef.ask(GetCurrentStatus).mapTo[StatusCheckResponse], scaled(20 seconds))
      assert(resp.ok == ok)

      val actual = resp.systems.toList.map(NameAndStatus.tupled).sortBy(_.name)
      val expected = systems.map { case (m, s) => m.name -> s } map { NameAndStatus.tupled } sortBy(_.name)
      assert(actual.length == expected.length)

      actual.zip(expected) map { case (a, e) =>
        assert(a.name == e.name)
        assert(a.status.ok == e.status.ok)
        val isEmpty = a.status.messages.isEmpty && e.status.messages.isEmpty
        lazy val actualPrefixedByExpected = for {
          am <- a.status.messages.toList.flatten.headOption
          em <- e.status.messages.toList.flatten.headOption
        } yield am.startsWith(em)
        assert(isEmpty || actualPrefixedByExpected.contains(true), s"Instead, a.status.messages = ${a.status.messages.toList.flatten} and e.status.messages = ${e.status.messages.toList.flatten}")
      } head
    }
  }

  "HealthMonitor" should "start with unknown status for all subsystems" in {
    val hm = TestActorRef(new TestHealthMonitorActor() {
      override def subsystems: Set[MonitoredSubsystem] = Set(SuccessSubsystem)
      // The default implementation of initialize kicks off an actual status check which may or may not complete before
      // the `ask` below comes back. What this test wants to know is the initial statuses of monitored subsystems
      // before the checks run, so this override neuters the initialization.
      override def initialize(): Unit = ()
    })
    eventualStatus(hm, ok = false, (SuccessSubsystem, UnknownStatus))
  }

  it should "return ok status when all subsystems succeed" in {
    val hm = TestActorRef(new TestHealthMonitorActor {
      override def subsystems: Set[MonitoredSubsystem] = Set(SuccessSubsystem)
    })
    eventualStatus(hm, ok = true, (SuccessSubsystem, OkStatus))
  }

  it should "fail if any subsystems fail but correctly bin them" in {
    val hm = TestActorRef(new TestHealthMonitorActor {
      override def subsystems: Set[MonitoredSubsystem] = Set(SuccessSubsystem, FailureSubsystem)
      /*
      Unlike success-statuses, failure-statuses retry a number of times before storing their values.
      While the failures-status is retrying, the success-status ends up becoming "stale", and flips to "unknown".
      So, increase the timeout for marking the success-status as stale and set the failureRetryCount to 0
       */
      override lazy val staleThreshold = scaled(10.seconds)
      override lazy val failureRetryCount = 0
    })
    eventualStatus(hm, ok = false, (SuccessSubsystem, OkStatus), (FailureSubsystem, FailedStatus))
  }

  it should "properly invalidate formerly good stale cache entries" in {
    val hm = TestActorRef(new TestHealthMonitorActor {
      override def subsystems: Set[MonitoredSubsystem] = Set(SuccessSubsystem)
      // a neutralized sweep check will cause staleness
      override private[healthmonitor] def scheduleSweepCheck(subsystem: MonitoredSubsystem): Unit = ()
    })
    eventualStatus(hm, ok = true, (SuccessSubsystem, OkStatus))
    eventualStatus(hm, ok = false, (SuccessSubsystem, UnknownStatus))
  }

  it should "properly invalidate formerly bad stale cache entries" in {
    val hm = TestActorRef(new TestHealthMonitorActor {
      override def subsystems: Set[MonitoredSubsystem] = Set(FailureSubsystem)
      // a neutralized sweep check will cause staleness
      override private[healthmonitor] def scheduleSweepCheck(subsystem: MonitoredSubsystem): Unit = ()
    })

    eventualStatus(hm, ok = false, (FailureSubsystem, FailedStatus))
    eventualStatus(hm, ok = false, (FailureSubsystem, UnknownStatus))
  }

  it should "handle timed out futures" in {
    var first = true
    def timeOutCheck(): Future[SubsystemStatus] = {
      if (first) {
        first = false
        Future.successful(OkStatus)
      } else {
        Future {
          Thread.sleep(scaled(10.seconds).toMillis)
          OkStatus
        }
      }
    }

    val timeoutSubsystem = MonitoredSubsystem("Timeout", () => timeOutCheck())

    val hm = TestActorRef(new TestHealthMonitorActor {
      override def subsystems: Set[MonitoredSubsystem] = Set(timeoutSubsystem)
      override lazy val failureRetryCount = 0
    })

    eventualStatus(hm, ok = true, (timeoutSubsystem, OkStatus))
    eventualStatus(hm, ok = false, (timeoutSubsystem, TimedOutStatus))
  }

  it should "retry failures the appropriate amount of times before failing" in {
    // Start up an actor with a status check function that's ok the first time but subsequently fails, incrementing a
    // counter each time it is invoked.
    var subsystemStatusCheckCount = 0
    var failedCount = 0
    def okThenFailThenHangStatusChecker(): Future[SubsystemStatus] = {
      subsystemStatusCheckCount = subsystemStatusCheckCount + 1
      subsystemStatusCheckCount match {
        case 1 =>
          Future.successful(SubsystemStatus(ok = true, messages = None))
        case 2 | 3 | 4 | 5 =>
          failedCount = failedCount + 1
          Future.failed(new RuntimeException("womp womp"))
        case _ =>
          Future {
            Thread.sleep(scaled(5.seconds).toMillis)
            SubsystemStatus(ok = true, messages = None)
          }
      }
    }

    var statusStores = List.empty[SubsystemStatus]

    val failAndCountSubsystem = MonitoredSubsystem("FailAndCountSubsystem", () => okThenFailThenHangStatusChecker())

    val hm = TestActorRef(new TestHealthMonitorActor {
      override def subsystems: Set[MonitoredSubsystem] = Set(failAndCountSubsystem)

      override def store(monitoredSubsystem: MonitoredSubsystem, subsystemStatus: SubsystemStatus): Unit = {
        statusStores = statusStores :+ subsystemStatus
        super.store(monitoredSubsystem, subsystemStatus)
      }
    })

    eventualStatus(hm, ok = true, (failAndCountSubsystem, OkStatus))
    eventualStatus(hm, ok = false, (failAndCountSubsystem, FailedStatus))
    eventualStatus(hm, ok = false, (failAndCountSubsystem, TimedOutStatus))

    // The initial failed try plus 3 failed retries
    assert(failedCount == 4)
    // OK, failed with womp womp, failed with Timed out
    assert(statusStores.size == 3)
    val statusStoreMessages: List[Option[String]] = statusStores.map(_.messages.flatMap(_.headOption))

    val expectedMessagePrefixes = List(None, Option("womp womp"), Option("Timed out"))
    statusStoreMessages.zip(expectedMessagePrefixes) foreach {
      case (a, e) => assert(a.isEmpty && e.isEmpty || a.exists(_.startsWith(e.get)))
    }
  }
}

object ProtoHealthMonitorServiceActorSpec {
  val SuccessSubsystem = MonitoredSubsystem("Success", () => mockCheckSuccess())
  val FailureSubsystem = MonitoredSubsystem("Failure", () => mockCheckFailure())
  val MultipleSubsystems = Set(SuccessSubsystem, FailureSubsystem)

  val FailedStatus = SubsystemStatus(ok = false, Option(List("womp womp")))
  val TimedOutStatus = SubsystemStatus(ok = false, Option(List("Timed out")))

  abstract class TestHealthMonitorActor(override val serviceConfig: Config = ConfigFactory.empty())
    extends ProtoHealthMonitorServiceActor with ScaledTimeSpans {
    override lazy val staleThreshold = scaled(3.seconds)
    override lazy val failureRetryInterval = scaled(100.milliseconds)
    override lazy val sweepInterval = scaled(200.milliseconds)
    override lazy val futureTimeout = scaled(1.second)
  }

  private def mockCheckSuccess(): Future[SubsystemStatus] = {
    Future.successful(OkStatus)
  }

  private def mockCheckFailure(): Future[SubsystemStatus] = {
    // Pay no mind, just needed to make sure the checks for both subsystems in the "binning" test get to run.
    Thread.`yield`()
    Future.failed(new RuntimeException("womp womp"))
  }
}
