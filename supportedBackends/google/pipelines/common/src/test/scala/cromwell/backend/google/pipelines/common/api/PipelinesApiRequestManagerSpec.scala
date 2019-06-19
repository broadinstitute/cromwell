package cromwell.backend.google.pipelines.common.api

import akka.actor.{ActorLogging, ActorRef, Props}
import akka.testkit.{TestActorRef, TestProbe, _}
import cats.data.NonEmptyList
import com.google.api.client.googleapis.batch.BatchRequest
import cromwell.backend.BackendSingletonActorAbortWorkflow
import cromwell.backend.google.pipelines.common.PipelinesApiTestConfig
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager._
import cromwell.backend.google.pipelines.common.api.TestPipelinesApiRequestManagerSpec._
import cromwell.backend.standard.StandardAsyncJob
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.util.AkkaTestUtil
import cromwell.util.AkkaTestUtil.DeathTestActor
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric._
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpecLike, Matchers}

import scala.collection.immutable.Queue
import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.util.Random

class PipelinesApiRequestManagerSpec extends TestKitSuite("PipelinesApiRequestManagerSpec") with FlatSpecLike with Matchers with Eventually {

  behavior of "PipelinesApiRequestManager"

  implicit val TestExecutionTimeout = 10.seconds.dilated
  implicit val DefaultPatienceConfig = PatienceConfig(TestExecutionTimeout)
  val AwaitAlmostNothing = 30.milliseconds.dilated
  val BatchSize = 5
  val registryProbe = TestProbe().ref
  val workflowId = WorkflowId.randomId()

  private def makePollRequest(snd: ActorRef, jobId: StandardAsyncJob) = new PAPIStatusPollRequest(workflowId, snd, null, jobId) {
    override def contentLength = 0
  }

  private def makeCreateRequest(contentSize: Long, snd: ActorRef, workflowId: WorkflowId = workflowId) = new PAPIRunCreationRequest(workflowId, snd, null) {
    override def contentLength = contentSize
  }

  it should "queue up and dispense status poll requests, in order" in {
    val statusPoller = TestProbe(name = "StatusPoller")
    val jaqmActor: TestActorRef[TestPipelinesApiRequestManager] = TestActorRef(TestPipelinesApiRequestManager.props(registryProbe, statusPoller.ref))

    var statusRequesters = ((0 until BatchSize * 2) map { i => i -> TestProbe(name = s"StatusRequester_$i") }).toMap

    // Initially, we should have no work:
    jaqmActor.tell(msg = PipelinesApiRequestManager.PipelinesWorkerRequestWork(BatchSize), sender = statusPoller.ref)
    statusPoller.expectMsg(max = TestExecutionTimeout, obj = PipelinesApiRequestManager.NoWorkToDo)

    // Send a few status poll requests:
    statusRequesters foreach { case (index, probe) =>
      jaqmActor.tell(msg = makePollRequest(probe.ref, StandardAsyncJob(index.toString)), sender = probe.ref)
    }

    // Should have no messages to the actual statusPoller yet:
    statusPoller.expectNoMessage(max = AwaitAlmostNothing)

    // Verify batches:
    2 times {
      jaqmActor.tell(msg = PipelinesApiRequestManager.PipelinesWorkerRequestWork(BatchSize), sender = statusPoller.ref)
      statusPoller.expectMsgPF(max = TestExecutionTimeout) {
        case PipelinesApiRequestManager.PipelinesApiWorkBatch(workBatch) =>
          val requesters = statusRequesters.take(BatchSize)
          statusRequesters = statusRequesters.drop(BatchSize)

          val zippedWithRequesters = workBatch.toList.zip(requesters)
          zippedWithRequesters foreach { case (pollQuery, (index, testProbe)) =>
            pollQuery.requester should be(testProbe.ref)
            pollQuery.asInstanceOf[PAPIStatusPollRequest].jobId should be(StandardAsyncJob(index.toString))
          }
        case other => fail(s"Unexpected message: $other")
      }
    }

    // Finally, we should have no work:
    jaqmActor.tell(msg = PipelinesApiRequestManager.PipelinesWorkerRequestWork(BatchSize), sender = statusPoller.ref)
    statusPoller.expectMsg(max = TestExecutionTimeout, obj = PipelinesApiRequestManager.NoWorkToDo)

    jaqmActor.underlyingActor.testPollerCreations should be(1)
  }

  it should "reject create requests above maxBatchSize" in {
    val statusPoller = TestProbe(name = "StatusPoller")
    val jaqmActor: TestActorRef[TestPipelinesApiRequestManager] = TestActorRef(TestPipelinesApiRequestManager.props(registryProbe, statusPoller.ref))

    val statusRequester = TestProbe()

    // Send a create request
    val request = makeCreateRequest(15 * 1024 * 1024, statusRequester.ref)
    jaqmActor.tell(msg = request, sender = statusRequester.ref)

    statusRequester.expectMsgClass(classOf[PipelinesApiRunCreationQueryFailed])

    jaqmActor.underlyingActor.queueSize shouldBe 0
  }

  it should "respect the maxBatchSize when beheading the queue" in {
    val statusPoller = TestProbe(name = "StatusPoller")
    // maxBatchSize is 14MB, which mean we can take 2 queries of 5MB but not 3
    val jaqmActor: TestActorRef[TestPipelinesApiRequestManager] = TestActorRef(TestPipelinesApiRequestManager.props(registryProbe, statusPoller.ref))

    val statusRequester = TestProbe()

    // Enqueue 3 create requests
    1 to 3 foreach { _ =>
      val request = makeCreateRequest(5 * 1024 * 1024, statusRequester.ref)
      jaqmActor.tell(msg = request, sender = statusRequester.ref)
    }

    // ask for a batch
    jaqmActor.tell(msg = PipelinesApiRequestManager.PipelinesWorkerRequestWork(BatchSize), sender = statusPoller.ref)

    // We should get only 2 requests back
    statusPoller.expectMsgPF(max = TestExecutionTimeout) {
      case PipelinesApiRequestManager.PipelinesApiWorkBatch(workBatch) => workBatch.toList.size shouldBe 2
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

      val statusPoller1 = TestActorRef[TestPapiWorkerActor](Props(new TestPapiWorkerActor(1)), TestActorRef(new AkkaTestUtil.StoppingSupervisor()), "statusPoller1")
      val statusPoller2 = TestActorRef[TestPapiWorkerActor](Props(new TestPapiWorkerActor(1)), TestActorRef(new AkkaTestUtil.StoppingSupervisor()), "statusPoller2")
      val statusPoller3 = TestActorRef[TestPapiWorkerActor](Props(new TestPapiWorkerActor(1)), TestActorRef(new AkkaTestUtil.StoppingSupervisor()), "statusPoller3")
      val jaqmActor: TestActorRef[TestPipelinesApiRequestManager] = TestActorRef(TestPipelinesApiRequestManager.props(registryProbe, statusPoller1, statusPoller2, statusPoller3), s"TestJesApiQueryManager-${Random.nextInt()}")

      val emptyActor = system.actorOf(Props.empty)

      // Send a few status poll requests:
      BatchSize indexedTimes { index =>
        val request = makePollRequest(emptyActor, StandardAsyncJob(index.toString))
        jaqmActor.tell(msg = request, sender = emptyActor)
      }

      // The work queue should be filled and the manager should have one poller active:
      eventually {
        jaqmActor.underlyingActor.queueSize should be (BatchSize)
        jaqmActor.underlyingActor.testPollerCreations should be (1)
        jaqmActor.underlyingActor.statusPollers.size should be(1)
        jaqmActor.underlyingActor.statusPollers.head should be(statusPoller1)
      }

      // Status poller 1 should get some work:
      jaqmActor ! SpecOnly_TriggerRequestForWork(jaqmActor, BatchSize)

      // Therefore, there's no work left in the work queue because it's all been given to the worker:
      eventually {
        jaqmActor.underlyingActor.queueSize should be (0)
        statusPoller1.underlyingActor.workToDo.map(_.size) should be(Some(BatchSize))
      }

      // Uh oh, something happens to status poller 1:
      stopMethod(statusPoller1)

      // The work queue receives all the work that couldn't be completed, and a new poller is created:
      eventually {
        jaqmActor.underlyingActor.queueSize should be (BatchSize)
        jaqmActor.underlyingActor.testPollerCreations should be (2)
        jaqmActor.underlyingActor.statusPollers.size should be(1)
        jaqmActor.underlyingActor.statusPollers.head should be(statusPoller2)
      }

      // Next, status poller 2 should be able to get some work:
      jaqmActor ! SpecOnly_TriggerRequestForWork(jaqmActor, BatchSize)

      // The queue is emptied again and this time statusPoller2 has it:
      eventually {
        jaqmActor.underlyingActor.queueSize should be (0)
        statusPoller2.underlyingActor.workToDo.map(_.size) should be(Some(BatchSize))
      }

      // Uh oh, something *also* happens to status poller 2!!
      stopMethod(statusPoller2)

      // The work queue receives all the work that couldn't be completed *again*, and a new poller is created *again*:
      eventually {
        jaqmActor.underlyingActor.queueSize should be (BatchSize)
        jaqmActor.underlyingActor.testPollerCreations should be (3)
        jaqmActor.underlyingActor.statusPollers.size should be(1)
        jaqmActor.underlyingActor.statusPollers.head should be(statusPoller3)
      }
    }
  }

  it should "remove run requests from queue when receiving an abort message" in {
    val statusPoller = TestProbe(name = "StatusPoller")
    // maxBatchSize is 14MB, which mean we can take 2 queries of 5MB but not 3
    val jaqmActor: TestActorRef[TestPipelinesApiRequestManager] = TestActorRef(TestPipelinesApiRequestManager.props(registryProbe, statusPoller.ref))

    // Enqueue 3 create requests
    val workflowIdA = WorkflowId.randomId()
    val workflowIdB = WorkflowId.randomId()
    jaqmActor ! makeCreateRequest(0, emptyActor, workflowIdA)
    jaqmActor ! makeCreateRequest(0, emptyActor, workflowIdA)
    jaqmActor ! makeCreateRequest(0, emptyActor, workflowIdB)

    // abort workflow A
    jaqmActor ! BackendSingletonActorAbortWorkflow(workflowIdA)

    // It should remove all and only run requests for workflow A
    eventually {
      jaqmActor.underlyingActor.queueSize shouldBe 1
      jaqmActor.underlyingActor.workQueue.head.asInstanceOf[PipelinesApiRequestManager.PAPIRunCreationRequest].workflowId shouldBe workflowIdB
    }
  }
}

object TestPipelinesApiRequestManagerSpec {
  implicit class intWithTimes(n: Int) {
    def times(f: => Unit) = 1 to n foreach { _ => f }
    def indexedTimes(f: Int => Unit) = 0 until n foreach { i => f(i) }
  }
}

/**
  * This test class allows us to hook into the JesApiQueryManager's makeStatusPoller and provide our own TestProbes instead
  */
class TestPipelinesApiRequestManager(qps: Int Refined Positive,
                                     requestWorkers: Int Refined Positive,
                                     registry: ActorRef,
                                     availableRequestWorkers: ActorRef*)
  extends PipelinesApiRequestManager(qps, requestWorkers, registry)(new MockPipelinesRequestHandler) {

  var testProbeQueue: Queue[ActorRef] = _
  var testPollerCreations: Int = _

  private def askForWorkReceive: PartialFunction[Any, Unit] = {
    case afw: SpecOnly_TriggerRequestForWork => statusPollers.head ! afw
  }

  override def receive: PartialFunction[Any, Unit] = askForWorkReceive orElse super.receive

  private def init() = {
    testProbeQueue = Queue(availableRequestWorkers: _*)
    testPollerCreations = 0
  }

  override private[api] lazy val nbWorkers = 1
  override private[api] def resetAllWorkers() = {
    val pollers = Vector.fill(1) { makeAndWatchWorkerActor() }
    pollers
  }

  override private[api] def makeWorkerActor(): ActorRef = {
    // Initialize the queue, if necessary:
    if (testProbeQueue == null) {
      init()
    }

    // Register that the creation was requested:
    testPollerCreations += 1

    // Pop the queue to get the next test probe:
    val (probe, newQueue) = testProbeQueue.dequeue
    testProbeQueue = newQueue
    probe
  }

  def queueSize = workQueue.size
  def statusPollerEquals(otherStatusPoller: ActorRef) = statusPollers sameElements Array(otherStatusPoller)
}

/**
  * Grabs some work but doesn't ever work on it (let alone completing it!)
  */
class TestPapiWorkerActor(maxBatchSize: Int) extends DeathTestActor with ActorLogging {

  var workToDo: Option[NonEmptyList[PAPIApiRequest]] = None

  override def receive = stoppingReceive orElse {

    case PipelinesApiWorkBatch(work) =>
      log.info(s"received ${work.size} requests from ${sender().path.name}")
      workToDo = Some(work)
    case NoWorkToDo =>
      log.info(s"received 0 requests from ${sender().path.name}")
      workToDo = None
    case SpecOnly_TriggerRequestForWork(manager, batchSize) => requestWork(manager, batchSize)
  }

  def requestWork(papiRequestManager: ActorRef, batchSize: Int): Unit = {
    papiRequestManager ! PipelinesWorkerRequestWork(batchSize)
  }
}

case class SpecOnly_TriggerRequestForWork(papiRequestManager: ActorRef, batchSize: Int)

class MockPipelinesRequestHandler extends PipelinesApiRequestHandler {
  override def makeBatchRequest = throw new UnsupportedOperationException
  override def enqueue[T <: PipelinesApiRequestManager.PAPIApiRequest](papiApiRequest: T, batchRequest: BatchRequest, pollingManager: ActorRef)
                                                                      (implicit ec: ExecutionContext)= throw new UnsupportedOperationException
}

object TestPipelinesApiRequestManager {
  import PipelinesApiTestConfig._

  def props(registryProbe: ActorRef, statusPollers: ActorRef*): Props = {
    Props(new TestPipelinesApiRequestManager(
      papiConfiguration.papiAttributes.qps,
      papiConfiguration.papiAttributes.requestWorkers,
      registryProbe,
      statusPollers: _*)
    )
  }
}
