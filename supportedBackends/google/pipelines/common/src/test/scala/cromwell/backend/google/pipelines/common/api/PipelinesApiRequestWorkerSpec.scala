package cromwell.backend.google.pipelines.common.api

import akka.actor.{ActorRef, Props}
import akka.testkit.{TestActorRef, TestProbe, _}
import cats.data.NonEmptyList
import com.google.api.client.googleapis.batch.BatchRequest
import com.google.api.client.googleapis.batch.json.JsonBatchCallback
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpRequest
import cromwell.backend.google.pipelines.common.PipelinesApiConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.{PAPIApiException, PAPIApiRequestFailed, PAPIStatusPollRequest, PipelinesWorkerRequestWork}
import cromwell.backend.google.pipelines.common.api.TestPipelinesApiRequestWorker.{CallbackFailure, CallbackSuccess, PipelinesApiBatchCallbackResponse}
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

abstract class PipelinesApiRequestWorkerSpec[O >: Null]
  extends TestKitSuite("PipelinesApiRequestWorker") with FlatSpecLike with Matchers with Eventually with BeforeAndAfter with Mockito {

  implicit var batchHandler: TestPipelinesApiBatchHandler[O]
  
  behavior of "PipelinesApiRequestWorker"

  implicit val TestExecutionTimeout = 10.seconds.dilated
  implicit val DefaultPatienceConfig = PatienceConfig(TestExecutionTimeout)
  val AwaitAlmostNothing = 30.milliseconds.dilated

  var managerProbe: TestProbe = _
  var workerActor: TestActorRef[TestPipelinesApiRequestWorker] = _

  it should "correctly calculate batch intervals" in {
    import eu.timepit.refined.auto._
    PipelinesApiRequestManager.determineBatchInterval(10) should be(11111.milliseconds)
    PipelinesApiRequestManager.determineBatchInterval(100000) shouldBe 1.millisecond
  }

  it should "query for work and wait for a reply" in {
    managerProbe.expectMsgClass(max = TestExecutionTimeout, c = classOf[PipelinesApiRequestManager.PipelinesWorkerRequestWork])
    managerProbe.expectNoMessage(max = AwaitAlmostNothing)
  }

  it should "respond correctly with various run statuses" in {
    managerProbe.expectMsgClass(max = TestExecutionTimeout, c = classOf[PipelinesApiRequestManager.PipelinesWorkerRequestWork])

    val requester1 = TestProbe()
    val query1 = PAPIStatusPollRequest(null, requester1.ref, null, null)
    val requester2 = TestProbe()
    val query2 = PAPIStatusPollRequest(null, requester2.ref, null, null)
    val requester3 = TestProbe()
    val query3 = PAPIStatusPollRequest(null, requester3.ref, null, null)

    // For two requests the callback succeeds (first with RunStatus.Success, then RunStatus.Failed). The third callback fails (simulating a network timeout, for example):
    batchHandler.callbackResponses :+= CallbackSuccess
    batchHandler.callbackResponses :+= CallbackSuccess
    batchHandler.callbackResponses :+= CallbackFailure

    val successStatus = RunStatus.Success(Seq.empty[ExecutionEvent], None, None, None)
    val failureStatus = RunStatus.UnsuccessfulRunStatus(Status.UNIMPLEMENTED, Option.empty[String], Seq.empty[ExecutionEvent], None, None, None, false)
    batchHandler.operationStatusResponses :+= successStatus
    batchHandler.operationStatusResponses :+= failureStatus

    workerActor.tell(msg = PipelinesApiRequestManager.PipelinesApiWorkBatch(NonEmptyList(query1, List(query2, query3))), sender = managerProbe.ref)
    eventually { batchHandler.runBatchRequested should be(true) }

    // The manager shouldn't have been asked for more work yet:
    managerProbe.expectNoMessage(max = AwaitAlmostNothing)

    // Ok, let's trigger the callbacks:
    batchHandler.executeBatch()

    requester1.expectMsg(successStatus)
    requester2.expectMsg(failureStatus)
    requester3.expectNoMessage(max = AwaitAlmostNothing)

    // Requester3 expected nothing... Instead, the manager expects an API failure notification and then a request for more work:
    managerProbe.expectMsgPF(TestExecutionTimeout) {
      case failure: PAPIApiRequestFailed =>
        if (!failure.cause.isInstanceOf[PAPIApiException]) fail("Unexpected failure cause class: " + failure.cause.getClass.getSimpleName)
        if (failure.query != query2 && failure.query != query3) fail("Unexpected query caused failure: " + failure.query)
    }
    managerProbe.expectMsg(PipelinesWorkerRequestWork(PipelinesApiRequestWorker.MaxBatchSize))
    managerProbe.expectNoMessage(max = AwaitAlmostNothing)
  }
}

class TestPipelinesApiRequestWorker(manager: ActorRef, qps: Int Refined Positive, registryProbe: ActorRef)(implicit batchHandler: TestPipelinesApiBatchHandler[_])
  extends PipelinesApiRequestWorker(manager, 10.milliseconds, registryProbe) with Mockito {
  override def createBatch() = null
  override def runBatch(batch: BatchRequest) = batchHandler.runBatch()
}

abstract class TestPipelinesApiBatchHandler[O >: Null] extends PipelinesApiRequestHandler {
  var operationStatusResponses: Queue[RunStatus] = Queue.empty
  var resultHandlers: Queue[JsonBatchCallback[O]] = Queue.empty
  var callbackResponses: Queue[PipelinesApiBatchCallbackResponse] = Queue.empty
  var runBatchRequested: Boolean = false

  def createBatch(genomicsInterface: O): BatchRequest = null
  def runBatch(): Unit = runBatchRequested = true

  def executeBatch(): Unit = {
    resultHandlers.zip(callbackResponses) foreach { case (handler, response) => response match {
      case CallbackSuccess => handler.onSuccess(null, null)
      case CallbackFailure =>
        val error: GoogleJsonError = new GoogleJsonError()
        handler.onFailure(error, null)
    }}
  }

  def enqueueStatusPollInBatch(pollingRequest: PAPIStatusPollRequest, batch: BatchRequest): Future[Try[Unit]] = {
    val completionPromise = Promise[Try[Unit]]()
    val resultHandler = statusPollResultHandler(pollingRequest, completionPromise)
    addStatusPollToBatch(null, batch, resultHandler)
    completionPromise.future
  }
  
  def statusPollResultHandler(pollRequest: PAPIStatusPollRequest, completionPromise: Promise[Try[Unit]]): JsonBatchCallback[O]

  def addStatusPollToBatch(httpRequest: HttpRequest, batch: BatchRequest, resultHandler: JsonBatchCallback[O]): Unit = resultHandlers :+= resultHandler
  
  def mockStatusInterpreter(operation: O): RunStatus = {
    val (status, newQueue) = operationStatusResponses.dequeue
    operationStatusResponses = newQueue
    status
  }
}

object TestPipelinesApiRequestWorker {
  def props(manager: ActorRef, jesConfiguration: PipelinesApiConfiguration, registryProbe: ActorRef)
           (implicit batchHandler: TestPipelinesApiBatchHandler[_]): Props = {
    Props(new TestPipelinesApiRequestWorker(
      manager,
      jesConfiguration.papiAttributes.qps,
      registryProbe
    ))
  }
  
  sealed trait PipelinesApiBatchCallbackResponse
  case object CallbackSuccess extends PipelinesApiBatchCallbackResponse
  case object CallbackFailure extends PipelinesApiBatchCallbackResponse
}
