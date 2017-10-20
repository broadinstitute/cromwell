package cromwell.backend.impl.jes.statuspolling

import akka.actor.{ActorRef, Props}
import akka.testkit.{TestActorRef, TestProbe, _}
import cats.data.NonEmptyList
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpRequest
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.Operation
import cromwell.backend.impl.jes.statuspolling.JesApiQueryManager.{JesApiException, JesApiQueryFailed, JesStatusPollQuery, RequestJesPollingWork}
import cromwell.backend.impl.jes.statuspolling.TestJesPollingActor.{CallbackFailure, CallbackSuccess, JesBatchCallbackResponse}
import cromwell.backend.impl.jes.{JesConfiguration, Run, RunStatus}
import cromwell.core.{ExecutionEvent, TestKitSuite}
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Positive
import io.grpc.Status
import org.scalatest.concurrent.Eventually
import org.scalatest.{BeforeAndAfter, FlatSpecLike, Matchers}
import org.specs2.mock.Mockito

import scala.collection.immutable.Queue
import scala.concurrent.duration._
import scala.concurrent.{Future, Promise}
import scala.util.Try

class JesPollingActorSpec extends TestKitSuite("JesPollingActor") with FlatSpecLike with Matchers with Eventually with BeforeAndAfter with Mockito {

  behavior of "JesPollingActor"

  implicit val TestExecutionTimeout = 10.seconds.dilated
  implicit val DefaultPatienceConfig = PatienceConfig(TestExecutionTimeout)
  val AwaitAlmostNothing = 30.milliseconds.dilated

  import cromwell.backend.impl.jes.JesTestConfig.JesBackendConfigurationDescriptor
  val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

  var managerProbe: TestProbe = _
  var jpActor: TestActorRef[TestJesPollingActor] = _

  it should "correctly calculate batch intervals" in {
    import eu.timepit.refined.auto._
    JesPollingActor.determineBatchInterval(10) should be(11.seconds)
    JesPollingActor.determineBatchInterval(100000) shouldBe 1.seconds
  }

  it should "query for work and wait for a reply" in {
    managerProbe.expectMsgClass(max = TestExecutionTimeout, c = classOf[JesApiQueryManager.RequestJesPollingWork])
    managerProbe.expectNoMsg(max = AwaitAlmostNothing)
  }

  it should "respond correctly with various run statuses" in {
    managerProbe.expectMsgClass(max = TestExecutionTimeout, c = classOf[JesApiQueryManager.RequestJesPollingWork])

    val requester1 = TestProbe()
    val query1 = JesStatusPollQuery(requester1.ref, Run(null, null))
    val requester2 = TestProbe()
    val query2 = JesStatusPollQuery(requester2.ref, Run(null, null))
    val requester3 = TestProbe()
    val query3 = JesStatusPollQuery(requester3.ref, Run(null, null))

    // For two requests the callback succeeds (first with RunStatus.Success, then RunStatus.Failed). The third callback fails (simulating a network timeout, for example):
    jpActor.underlyingActor.callbackResponses :+= CallbackSuccess
    jpActor.underlyingActor.callbackResponses :+= CallbackSuccess
    jpActor.underlyingActor.callbackResponses :+= CallbackFailure

    val successStatus = RunStatus.Success(Seq.empty[ExecutionEvent], None, None, None)
    val failureStatus = RunStatus.UnsuccessfulRunStatus(Status.UNIMPLEMENTED, Option.empty[String], Seq.empty[ExecutionEvent], None, None, None)
    jpActor.underlyingActor.operationStatusResponses :+= successStatus
    jpActor.underlyingActor.operationStatusResponses :+= failureStatus

    jpActor.tell(msg = JesApiQueryManager.JesPollingWorkBatch(NonEmptyList(query1, List(query2, query3))), sender = managerProbe.ref)
    eventually { jpActor.underlyingActor.resultHandlers.size should be(3) }
    eventually { jpActor.underlyingActor.runBatchRequested should be(true) }

    // The manager shouldn't have been asked for more work yet:
    managerProbe.expectNoMsg(max = AwaitAlmostNothing)

    // Ok, let's trigger the callbacks:
    jpActor.underlyingActor.executeBatch()

    requester1.expectMsg(successStatus)
    requester2.expectMsg(failureStatus)
    requester3.expectNoMsg(max = AwaitAlmostNothing)

    // Requester3 expected nothing... Instead, the manager expects an API failure notification and then a request for more work:
    managerProbe.expectMsgPF(TestExecutionTimeout) {
      case failure: JesApiQueryFailed =>
        if (!failure.cause.isInstanceOf[JesApiException]) fail("Unexpected failure cause class: " + failure.cause.getClass.getSimpleName)
        if (failure.query != query2 && failure.query != query3) fail("Unexpected query caused failure: " + failure.query)
    }
    managerProbe.expectMsg(RequestJesPollingWork(JesPollingActor.MaxBatchSize))
    managerProbe.expectNoMsg(max = AwaitAlmostNothing)
  }

  before {
    managerProbe = TestProbe()
    val registryProbe = TestProbe()
    jpActor = TestActorRef(TestJesPollingActor.props(managerProbe.ref, jesConfiguration, registryProbe.ref), managerProbe.ref)
  }
}

/**
  * Testable JES polling actor.
  * - Mocks out the methods which actually call out to JES, and allows the callbacks to be triggered in a testable way
  * - Also waits a **lot** less time before polls!
  */
class TestJesPollingActor(manager: ActorRef, qps: Int Refined Positive, registryProbe: ActorRef) extends JesPollingActor(manager, qps, registryProbe) with Mockito {

  override lazy val batchInterval = 10.milliseconds

  var operationStatusResponses: Queue[RunStatus] = Queue.empty
  var resultHandlers: Queue[JsonBatchCallback[Operation]] = Queue.empty
  var callbackResponses: Queue[JesBatchCallbackResponse] = Queue.empty
  var runBatchRequested: Boolean = false

  override private[statuspolling] def createBatch(genomicsInterface: Genomics): BatchRequest = null
  override private[statuspolling] def runBatch(batch: BatchRequest): Unit = runBatchRequested = true

  def executeBatch(): Unit = {
    resultHandlers.zip(callbackResponses) foreach { case (handler, response) => response match {
      case CallbackSuccess => handler.onSuccess(null, null)
      case CallbackFailure =>
        val error: GoogleJsonError = new GoogleJsonError()
        handler.onFailure(error, null)
    }}
  }

  override private [statuspolling] def enqueueStatusPollInBatch(pollingRequest: JesStatusPollQuery, batch: BatchRequest): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = statusPollResultHandler(pollingRequest, completionPromise)
    addStatusPollToBatch(null, batch, resultHandler)
    completionPromise.future
  }
  
  override def addStatusPollToBatch(httpRequest: HttpRequest, batch: BatchRequest, resultHandler: JsonBatchCallback[Operation]): Unit = resultHandlers :+= resultHandler
  override def interpretOperationStatus(operation: Operation): RunStatus = {
    val (status, newQueue) = operationStatusResponses.dequeue
    operationStatusResponses = newQueue
    status
  }
  override def mkErrorString(e: GoogleJsonError) = "NA"
}

object TestJesPollingActor {
  def props(manager: ActorRef, jesConfiguration: JesConfiguration, registryProbe: ActorRef) = Props(new TestJesPollingActor(manager, jesConfiguration.qps, registryProbe))

  sealed trait JesBatchCallbackResponse
  case object CallbackSuccess extends JesBatchCallbackResponse
  case object CallbackFailure extends JesBatchCallbackResponse
}
