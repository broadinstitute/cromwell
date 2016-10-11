package cromwell.backend.impl.jes.statuspolling

import akka.actor.{ActorRef, Props}
import akka.testkit.{TestActorRef, TestProbe}
import cromwell.backend.impl.jes.Run
import cromwell.core.TestKitSuite
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import akka.testkit._
import JesApiQueryManagerSpec._
import cromwell.util.AkkaTestUtil
import org.scalatest.concurrent.Eventually

import scala.collection.immutable.Queue

class JesApiQueryManagerSpec extends TestKitSuite("JesApiQueryManagerSpec") with FlatSpecLike with Matchers with Eventually {

  behavior of "JesApiQueryManagerSpec"

  implicit val TestExecutionTimeout = 10.seconds.dilated
  implicit val DefaultPatienceConfig = PatienceConfig(TestExecutionTimeout)
  val AwaitAlmostNothing = 30.milliseconds.dilated
  val BatchSize = 5

  it should "queue up and dispense status poll requests, in order" in {
    val statusPoller = TestProbe(name = "StatusPoller")
    val jaqmActor: TestActorRef[TestJesApiQueryManager] = TestActorRef(TestJesApiQueryManager.props(statusPoller.ref))

    var statusRequesters = ((0 until BatchSize * 2) map { i => i -> TestProbe(name = s"StatusRequester_$i") }).toMap

    // Initially, we should have no work:
    jaqmActor.tell(msg = JesApiQueryManager.RequestJesPollingWork(BatchSize), sender = statusPoller.ref)
    statusPoller.expectMsg(max = TestExecutionTimeout, obj = JesApiQueryManager.NoWorkToDo)

    // Send a few status poll requests:
    statusRequesters foreach { case (index, probe) =>
      jaqmActor.tell(msg = JesApiQueryManager.DoPoll(Run(index.toString, null)), sender = probe.ref)
    }

    // Should have no messages to the actual statusPoller yet:
    statusPoller.expectNoMsg(max = AwaitAlmostNothing)

    // Verify batches:
    2 times {
      jaqmActor.tell(msg = JesApiQueryManager.RequestJesPollingWork(BatchSize), sender = statusPoller.ref)
      statusPoller.expectMsgPF(max = TestExecutionTimeout) {
        case JesApiQueryManager.JesPollingWorkBatch(workBatch) =>
          val requesters = statusRequesters.take(BatchSize)
          statusRequesters = statusRequesters.drop(BatchSize)

          val zippedWithRequesters = workBatch.toList.zip(requesters)
          zippedWithRequesters foreach { case (pollQuery, (index, testProbe)) =>
            pollQuery.requester should be(testProbe.ref)
            pollQuery.run.runId should be(index.toString)
          }
      }
    }

    // Finally, we should have no work:
    jaqmActor.tell(msg = JesApiQueryManager.RequestJesPollingWork(BatchSize), sender = statusPoller.ref)
    statusPoller.expectMsg(max = TestExecutionTimeout, obj = JesApiQueryManager.NoWorkToDo)

    jaqmActor.underlyingActor.testPollerCreations should be(1)
  }

  AkkaTestUtil.actorDeathMethods(system) foreach { case (name, stopMethod) =>
    it should s"catch polling actors if they $name and then recreate them" in {

      val statusPoller1 = TestActorRef(Props(new AkkaTestUtil.DeathTestActor()), TestActorRef(new AkkaTestUtil.StoppingSupervisor()))
      val statusPoller2 = TestActorRef(Props(new AkkaTestUtil.DeathTestActor()))
      val jaqmActor: TestActorRef[TestJesApiQueryManager] = TestActorRef(TestJesApiQueryManager.props(statusPoller1, statusPoller2))

      val statusRequesters = ((0 until BatchSize * 2) map { i => i -> TestProbe(name = s"StatusRequester_$i") }).toMap

      // Send a few status poll requests:
      BatchSize indexedTimes { index =>
        val probe = statusRequesters(index)
        jaqmActor.tell(msg = JesApiQueryManager.DoPoll(Run(index.toString, null)), sender = probe.ref)
      }
      BatchSize indexedTimes { i =>
        val index = i + BatchSize // For the second half of the statusRequester set
        val probe = statusRequesters(index)
        jaqmActor.tell(msg = JesApiQueryManager.DoPoll(Run(index.toString, null)), sender = probe.ref)
      }

      // Request a set of work from the middle of the queue:
      val batchOffset = 2
      jaqmActor.tell(msg = JesApiQueryManager.RequestJesPollingWork(batchOffset), sender = statusPoller1)
      jaqmActor.tell(msg = JesApiQueryManager.RequestJesPollingWork(BatchSize), sender = statusPoller1)

      // Kill the original status poller:
      stopMethod(statusPoller1)

      // Only the appropriate requesters get an error:
      (0 until batchOffset) foreach { index =>
        val probe = statusRequesters(index)
        probe.expectNoMsg(max = AwaitAlmostNothing)
      }
      (batchOffset until batchOffset + BatchSize) foreach { index =>
        val probe = statusRequesters(index)
        probe.expectMsg(max = TestExecutionTimeout, hint = s"Polling error to requester #$index", obj = JesPollingActor.JesPollError)
      }
      (batchOffset + BatchSize until 2 * BatchSize) foreach { index =>
        val probe = statusRequesters(index)
        probe.expectNoMsg(max = AwaitAlmostNothing)
      }

      // Check the next status poller gets created:
      eventually { jaqmActor.underlyingActor.testPollerCreations should be(2) }
    }
  }
}

object JesApiQueryManagerSpec {
  implicit class intWithTimes(n: Int) {
    def times(f: => Unit) = 1 to n foreach { _ => f }
    def indexedTimes(f: Int => Unit) = 0 until n foreach { i => f(i) }
  }
}

/**
  * This test class allows us to hook into the JesApiQueryManager's makeStatusPoller and provide our own TestProbes instead
  */
class TestJesApiQueryManager(statusPollerProbes: ActorRef*) extends JesApiQueryManager {

  var testProbes: Queue[ActorRef] = _
  var testPollerCreations: Int = _

  private def init() = {
    testProbes = Queue(statusPollerProbes: _*)
    testPollerCreations = 0
  }

  override private[statuspolling] def makeStatusPoller(): ActorRef = {
    // Initialise the queue, if necessary:
    if (testProbes == null) {
      init()
    }

    // Register that the creation was requested:
    testPollerCreations += 1

    // Pop the queue to get the next test probe:
    val (probe, newQueue) = testProbes.dequeue
    testProbes = newQueue
    probe
  }
}

object TestJesApiQueryManager {
  def props(statusPollers: ActorRef*): Props = Props(new TestJesApiQueryManager(statusPollers: _*))
}
