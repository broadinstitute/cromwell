package cromwell.backend.impl.jes.statuspolling

import akka.actor.{ActorRef, Props}
import akka.testkit.{TestActorRef, TestProbe, _}
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.RunPipelineRequest
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.{JesApiRunCreationQueryFailed, JesRunCreationQuery, JesStatusPollQuery}
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManagerSpec._
import cromwell.backend.impl.jes.{JesConfiguration, Run}
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.TestKitSuite
import cromwell.util.AkkaTestUtil
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric._
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpecLike, Matchers}

import scala.collection.immutable.Queue
import scala.concurrent.duration._
import scala.util.Random

class JesApiQueryManagerSpec extends TestKitSuite("JesApiQueryManagerSpec") with FlatSpecLike with Matchers with Eventually {

  behavior of "JesApiQueryManagerSpec"

  implicit val TestExecutionTimeout = 10.seconds.dilated
  implicit val DefaultPatienceConfig = PatienceConfig(TestExecutionTimeout)
  val AwaitAlmostNothing = 30.milliseconds.dilated
  val BatchSize = 5
  val registryProbe = TestProbe().ref

  it should "queue up and dispense status poll requests, in order" in {
    val statusPoller = TestProbe(name = "StatusPoller")
    val jaqmActor: TestActorRef[TestJesApiQueryManager] = TestActorRef(TestJesApiQueryManager.props(100, registryProbe, statusPoller.ref))

    var statusRequesters = ((0 until BatchSize * 2) map { i => i -> TestProbe(name = s"StatusRequester_$i") }).toMap

    // Initially, we should have no work:
    jaqmActor.tell(msg = JesApiQueryManager.RequestJesPollingWork(BatchSize), sender = statusPoller.ref)
    statusPoller.expectMsg(max = TestExecutionTimeout, obj = JesApiQueryManager.NoWorkToDo)

    // Send a few status poll requests:
    statusRequesters foreach { case (index, probe) =>
      jaqmActor.tell(msg = JesApiQueryManager.DoPoll(Run(StandardAsyncJob(index.toString), null)), sender = probe.ref)
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
            pollQuery.asInstanceOf[JesStatusPollQuery].run.job should be(StandardAsyncJob(index.toString))
          }
        case other => fail(s"Unexpected message: $other")
      }
    }

    // Finally, we should have no work:
    jaqmActor.tell(msg = JesApiQueryManager.RequestJesPollingWork(BatchSize), sender = statusPoller.ref)
    statusPoller.expectMsg(max = TestExecutionTimeout, obj = JesApiQueryManager.NoWorkToDo)

    jaqmActor.underlyingActor.testPollerCreations should be(1)
  }

  it should "reject create requests above maxBatchSize" in {
    val statusPoller = TestProbe(name = "StatusPoller")
    val jaqmActor: TestActorRef[TestJesApiQueryManager] = TestActorRef(TestJesApiQueryManager.props(15 * 1024 * 1024, registryProbe, statusPoller.ref))

    val statusRequester = TestProbe()
    
    // Send a create request
    jaqmActor.tell(msg = JesApiQueryManager.DoCreateRun(null, null), sender = statusRequester.ref)

    statusRequester.expectMsgClass(classOf[JesApiRunCreationQueryFailed])

    jaqmActor.underlyingActor.queueSize shouldBe 0
  }

  it should "respect the maxBatchSize when beheading the queue" in {
    val statusPoller = TestProbe(name = "StatusPoller")
    // maxBatchSize is 14MB, which mean we can take 2 queries of 5MB but not 3
    val jaqmActor: TestActorRef[TestJesApiQueryManager] = TestActorRef(TestJesApiQueryManager.props(5 * 1024 * 1024, registryProbe, statusPoller.ref))

    val statusRequester = TestProbe()

    // Enqueue 3 create requests
    1 to 3 foreach { _ => jaqmActor.tell(msg = JesApiQueryManager.DoCreateRun(null, null), sender = statusRequester.ref) }

    // ask for a batch
    jaqmActor.tell(msg = JesApiQueryManager.RequestJesPollingWork(BatchSize), sender = statusPoller.ref)
    
    // We should get only 2 requests back
    statusPoller.expectMsgPF(max = TestExecutionTimeout) {
      case JesApiQueryManager.JesPollingWorkBatch(workBatch) => workBatch.toList.size shouldBe 2
      case other => fail(s"Unexpected message: $other")
    }
    
    // There should be 1 left in the queue
    jaqmActor.underlyingActor.queueSize shouldBe 1
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
      val jaqmActor: TestActorRef[TestJesApiQueryManager] = TestActorRef(TestJesApiQueryManager.props(100, registryProbe, statusPoller1, statusPoller2), s"TestJesApiQueryManager-${Random.nextInt()}")

      val emptyActor = system.actorOf(Props.empty)

       // Send a few status poll requests:
      BatchSize indexedTimes { index =>
        jaqmActor.tell(msg = JesApiQueryManager.DoPoll(Run(StandardAsyncJob(index.toString), null)), sender = emptyActor)
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
class TestJesApiQueryManager(qps: Int Refined Positive, createRequestSize: Long, registry: ActorRef, statusPollerProbes: ActorRef*) extends JesApiQueryManager(qps, registry) {
  var testProbes: Queue[ActorRef] = _
  var testPollerCreations: Int = _

  private def init() = {
    testProbes = Queue(statusPollerProbes: _*)
    testPollerCreations = 0
  }

  override private[statuspolling] def makeCreateQuery(replyTo: ActorRef, genomics: Genomics, rpr: RunPipelineRequest) = {
    new JesRunCreationQuery(replyTo, genomics, rpr) {
      override def contentLength = createRequestSize
    }
  }

  override private[statuspolling] def makePollQuery(replyTo: ActorRef, run: Run) = {
    new JesStatusPollQuery(replyTo, run) {
      override def contentLength = 0
    }
  }

  override private[statuspolling] def makeWorkerActor(): ActorRef = {
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

  def props(createRequestSize: Long, registryProbe: ActorRef, statusPollers: ActorRef*): Props = Props(new TestJesApiQueryManager(jesConfiguration.qps, createRequestSize, registryProbe, statusPollers: _*))
}
