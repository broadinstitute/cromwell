package cromwell.backend.impl.jes.statuspolling

import akka.actor.{ActorRef, Props}
import akka.testkit.{TestActorRef, TestProbe}
import cromwell.backend.impl.jes.{JesConfiguration, Run}
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
    /*
      This test creates two statusPoller ActorRefs which are handed to the TestJesApiQueryManager. Work is added to that query
      manager and then the first statusPoller requests work and is subsequently killed. The expectation is that:

      - The work will return to the workQueue of the query manager
      - The query manager will have registered a new statusPoller
      - That statusPoller is the second ActorRef (and artifact of TestJesApiQueryManager)
     */
    it should s"catch polling actors if they $name, recreate them and add work back to the queue" in {
      val statusPoller1 = TestActorRef(Props(new AkkaTestUtil.DeathTestActor()), TestActorRef(new AkkaTestUtil.StoppingSupervisor()))
      val statusPoller2 = TestActorRef(Props(new AkkaTestUtil.DeathTestActor()), TestActorRef(new AkkaTestUtil.StoppingSupervisor()))
      val jaqmActor: TestActorRef[TestJesApiQueryManager] = TestActorRef(TestJesApiQueryManager.props(statusPoller1, statusPoller2))

      val emptyActor = system.actorOf(Props.empty)

       // Send a few status poll requests:
      BatchSize indexedTimes { index =>
        jaqmActor.tell(msg = JesApiQueryManager.DoPoll(Run(index.toString, null)), sender = emptyActor)
      }

      jaqmActor.tell(msg = JesApiQueryManager.RequestJesPollingWork(BatchSize), sender = statusPoller1)

      stopMethod(statusPoller1)

      eventually {
        jaqmActor.underlyingActor.testPollerCreations should be (2)
        jaqmActor.underlyingActor.queueSize should be (BatchSize)
        jaqmActor.underlyingActor.statusPollerEquals(statusPoller2) should be (true)
      }
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
class TestJesApiQueryManager(qps: Int, statusPollerProbes: ActorRef*) extends JesApiQueryManager(qps) {
  var testProbes: Queue[ActorRef] = _
  var testPollerCreations: Int = _

  private def init() = {
    testProbes = Queue(statusPollerProbes: _*)
    testPollerCreations = 0
  }

  override private[statuspolling] def makeStatusPoller(): ActorRef = {
    // Initialize the queue, if necessary:
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

  def queueSize = workQueue.size
  def statusPollerEquals(otherStatusPoller: ActorRef) = statusPoller == otherStatusPoller
}

object TestJesApiQueryManager {
  import cromwell.backend.impl.jes.JesTestConfig.JesBackendConfigurationDescriptor
  val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

  def props(statusPollers: ActorRef*): Props = Props(new TestJesApiQueryManager(jesConfiguration.qps, statusPollers: _*))
}
